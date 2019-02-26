% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:25
% EndTime: 2019-02-26 20:59:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (135->34), mult. (167->45), div. (0->0), fcn. (193->8), ass. (0->26)
t12 = sin(qJ(3));
t15 = cos(qJ(3));
t23 = r_i_i_C(2) + qJ(5) + pkin(8);
t20 = r_i_i_C(3) + qJ(6);
t25 = pkin(5) + r_i_i_C(1);
t14 = cos(qJ(4));
t6 = pkin(4) * t14 + pkin(3);
t9 = qJ(4) + pkin(9);
t7 = sin(t9);
t8 = cos(t9);
t26 = t20 * t7 + t25 * t8 + t6;
t29 = t23 * t12 + t26 * t15;
t27 = t23 * t15;
t11 = sin(qJ(4));
t24 = pkin(4) * t11;
t13 = sin(qJ(1));
t22 = t12 * t13;
t16 = cos(qJ(1));
t21 = t12 * t16;
t19 = pkin(1) + pkin(7) + t24;
t18 = t12 * t6 + qJ(2) - t27;
t4 = -t13 * t7 + t8 * t21;
t3 = t13 * t8 + t7 * t21;
t2 = t16 * t7 + t8 * t22;
t1 = -t16 * t8 + t7 * t22;
t5 = [-t19 * t13 + t18 * t16 + t20 * t3 + t25 * t4, t13, t29 * t13, t20 * t2 - t25 * t1 + (-t11 * t22 + t14 * t16) * pkin(4), -t13 * t15, t1; t20 * t1 + t18 * t13 + t19 * t16 + t25 * t2, -t16, -t29 * t16, -t20 * t4 + t25 * t3 + (t11 * t21 + t13 * t14) * pkin(4), t16 * t15, -t3; 0, 0, -t12 * t26 + t27 (t20 * t8 - t25 * t7 - t24) * t15, t12, t15 * t7;];
Ja_transl  = t5;
