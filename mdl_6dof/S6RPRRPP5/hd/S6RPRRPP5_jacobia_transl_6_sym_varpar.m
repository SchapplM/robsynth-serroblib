% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RPRRPP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:38
% EndTime: 2019-02-26 20:58:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (151->28), mult. (181->36), div. (0->0), fcn. (211->7), ass. (0->24)
t19 = pkin(8) - r_i_i_C(3) - qJ(6);
t10 = pkin(9) + qJ(3);
t8 = sin(t10);
t17 = t19 * t8;
t9 = cos(t10);
t27 = t17 + pkin(3) * t9 + cos(pkin(9)) * pkin(2) + pkin(1);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t20 = pkin(4) + pkin(5) + r_i_i_C(1);
t21 = r_i_i_C(2) + qJ(5);
t26 = t21 * t12 + t20 * t14 + pkin(3);
t13 = sin(qJ(1));
t25 = t13 * t12;
t24 = t13 * t14;
t15 = cos(qJ(1));
t23 = t14 * t15;
t22 = t15 * t12;
t16 = t19 * t9 - t26 * t8;
t11 = -pkin(7) - qJ(2);
t4 = t9 * t23 + t25;
t3 = t9 * t22 - t24;
t2 = t9 * t24 - t22;
t1 = t9 * t25 + t23;
t5 = [-t21 * t1 - t11 * t15 - t27 * t13 - t20 * t2, t13, t16 * t15, -t20 * t3 + t21 * t4, t3, -t15 * t8; -t13 * t11 + t27 * t15 + t20 * t4 + t21 * t3, -t15, t16 * t13, -t20 * t1 + t21 * t2, t1, -t13 * t8; 0, 0, t26 * t9 + t17 (-t20 * t12 + t21 * t14) * t8, t8 * t12, t9;];
Ja_transl  = t5;
