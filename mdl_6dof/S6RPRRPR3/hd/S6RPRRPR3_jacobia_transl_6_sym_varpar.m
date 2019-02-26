% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:01
% EndTime: 2019-02-26 21:02:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (231->42), mult. (293->64), div. (0->0), fcn. (356->10), ass. (0->30)
t20 = cos(qJ(3));
t17 = sin(qJ(3));
t29 = -r_i_i_C(3) - pkin(9) + pkin(8);
t26 = t29 * t17;
t34 = t20 * pkin(3) + pkin(2) + t26;
t16 = sin(qJ(4));
t19 = cos(qJ(4));
t18 = cos(qJ(6));
t15 = sin(qJ(6));
t32 = pkin(4) + pkin(5);
t25 = t15 * r_i_i_C(2) - t32;
t22 = t18 * r_i_i_C(1) - t25;
t27 = -t15 * r_i_i_C(1) - qJ(5);
t23 = t18 * r_i_i_C(2) - t27;
t33 = t23 * t16 + t22 * t19 + pkin(3);
t14 = qJ(1) + pkin(10);
t12 = sin(t14);
t13 = cos(t14);
t31 = t16 * t20;
t4 = t12 * t31 + t13 * t19;
t3 = t4 * t18;
t30 = t19 * t20;
t24 = (-(t15 * t19 - t16 * t18) * r_i_i_C(1) - (t15 * t16 + t18 * t19) * r_i_i_C(2)) * t17;
t21 = -t33 * t17 + t29 * t20;
t7 = t12 * t16 + t13 * t30;
t6 = -t12 * t19 + t13 * t31;
t5 = t12 * t30 - t13 * t16;
t2 = t6 * t15 + t7 * t18;
t1 = -t7 * t15 + t6 * t18;
t8 = [-sin(qJ(1)) * pkin(1) + t13 * pkin(7) - t3 * r_i_i_C(2) + t27 * t4 - t22 * t5 - t34 * t12, 0, t21 * t13, -t22 * t6 + t23 * t7, t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); cos(qJ(1)) * pkin(1) + t12 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * qJ(5) + t32 * t7 + t34 * t13, 0, t21 * t12, -t3 * r_i_i_C(1) + t23 * t5 + t25 * t4, t4 (-t5 * t15 + t3) * r_i_i_C(1) + (-t4 * t15 - t5 * t18) * r_i_i_C(2); 0, 1, t33 * t20 + t26 (qJ(5) * t19 - t32 * t16) * t17 - t24, t17 * t16, t24;];
Ja_transl  = t8;
