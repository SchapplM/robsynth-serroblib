% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:48
% EndTime: 2019-02-26 20:58:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->26), mult. (144->34), div. (0->0), fcn. (167->7), ass. (0->24)
t24 = pkin(8) + r_i_i_C(2);
t10 = pkin(9) + qJ(3);
t8 = sin(t10);
t18 = t24 * t8;
t9 = cos(t10);
t27 = t18 + pkin(3) * t9 + cos(pkin(9)) * pkin(2) + pkin(1);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t19 = r_i_i_C(3) + qJ(5);
t25 = pkin(4) + r_i_i_C(1);
t26 = t19 * t12 + t25 * t14 + pkin(3);
t13 = sin(qJ(1));
t23 = t13 * t12;
t22 = t13 * t14;
t15 = cos(qJ(1));
t21 = t14 * t15;
t20 = t15 * t12;
t16 = t24 * t9 - t26 * t8;
t11 = -pkin(7) - qJ(2);
t4 = t9 * t21 + t23;
t3 = t9 * t20 - t22;
t2 = t9 * t22 - t20;
t1 = t9 * t23 + t21;
t5 = [-t19 * t1 - t11 * t15 - t27 * t13 - t25 * t2, t13, t16 * t15, t19 * t4 - t25 * t3, t3, 0; -t13 * t11 + t27 * t15 + t19 * t3 + t25 * t4, -t15, t16 * t13, -t25 * t1 + t19 * t2, t1, 0; 0, 0, t26 * t9 + t18 (-t25 * t12 + t19 * t14) * t8, t8 * t12, 0;];
Ja_transl  = t5;
