% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:51
% EndTime: 2019-02-26 22:08:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (258->48), mult. (532->80), div. (0->0), fcn. (675->12), ass. (0->36)
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t21 = pkin(11) + qJ(6);
t19 = sin(t21);
t20 = cos(t21);
t35 = pkin(5) * sin(pkin(11)) + qJ(4);
t31 = t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + t35;
t37 = pkin(3) + r_i_i_C(3) + pkin(10) + qJ(5);
t44 = t31 * t25 + t37 * t28 + pkin(2);
t43 = pkin(9) + cos(pkin(11)) * pkin(5) + pkin(4);
t42 = cos(qJ(1));
t23 = sin(pkin(6));
t27 = sin(qJ(1));
t41 = t23 * t27;
t40 = t23 * t28;
t29 = cos(qJ(2));
t39 = t23 * t29;
t38 = cos(pkin(6));
t36 = t23 * t42;
t34 = t27 * t38;
t33 = t38 * t42;
t32 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t43;
t26 = sin(qJ(2));
t12 = t26 * t33 + t27 * t29;
t3 = t12 * t25 + t28 * t36;
t4 = t12 * t28 - t25 * t36;
t14 = -t26 * t34 + t29 * t42;
t13 = t26 * t42 + t29 * t34;
t11 = t27 * t26 - t29 * t33;
t10 = t38 * t25 + t26 * t40;
t9 = t23 * t26 * t25 - t38 * t28;
t8 = t14 * t28 + t25 * t41;
t7 = t14 * t25 - t27 * t40;
t2 = t13 * t20 + t7 * t19;
t1 = -t13 * t19 + t7 * t20;
t5 = [-t27 * pkin(1) - t12 * pkin(2) + pkin(8) * t36 - t32 * t11 - t31 * t3 - t37 * t4, -t13 * t44 + t32 * t14, t31 * t8 - t37 * t7, t7, t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t14 * pkin(2) + pkin(8) * t41 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t43 * t13 + t35 * t7 + t37 * t8, -t11 * t44 + t32 * t12, -t37 * t3 + t31 * t4, t3, t4 (-t11 * t19 + t3 * t20) * r_i_i_C(1) + (-t11 * t20 - t3 * t19) * r_i_i_C(2); 0 (t32 * t26 + t44 * t29) * t23, t31 * t10 - t37 * t9, t9, t10 (t19 * t39 + t9 * t20) * r_i_i_C(1) + (-t9 * t19 + t20 * t39) * r_i_i_C(2);];
Ja_transl  = t5;
