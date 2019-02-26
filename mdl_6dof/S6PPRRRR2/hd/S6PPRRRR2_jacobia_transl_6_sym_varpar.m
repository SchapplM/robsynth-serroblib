% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:15
% EndTime: 2019-02-26 19:43:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (448->55), mult. (1116->98), div. (0->0), fcn. (1467->16), ass. (0->52)
t34 = cos(pkin(12));
t58 = sin(pkin(13));
t59 = sin(pkin(12));
t49 = t59 * t58;
t61 = cos(pkin(13));
t56 = t34 * t61;
t63 = cos(pkin(6));
t43 = -t63 * t56 + t49;
t33 = sin(pkin(6));
t60 = sin(pkin(7));
t57 = t33 * t60;
t62 = cos(pkin(7));
t71 = t34 * t57 + t43 * t62;
t50 = t59 * t61;
t55 = t34 * t58;
t44 = t63 * t50 + t55;
t54 = t59 * t33;
t70 = t44 * t62 - t60 * t54;
t24 = t63 * t55 + t50;
t37 = sin(qJ(3));
t65 = cos(qJ(3));
t11 = t24 * t37 + t71 * t65;
t32 = qJ(5) + qJ(6);
t30 = sin(t32);
t31 = cos(t32);
t12 = t24 * t65 - t71 * t37;
t64 = t34 * t33;
t19 = t43 * t60 - t62 * t64;
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t8 = t12 * t39 + t19 * t36;
t69 = (t11 * t31 - t8 * t30) * r_i_i_C(1) + (-t11 * t30 - t8 * t31) * r_i_i_C(2);
t25 = -t63 * t49 + t56;
t14 = t25 * t65 - t70 * t37;
t20 = t44 * t60 + t62 * t54;
t10 = t14 * t39 + t20 * t36;
t13 = t25 * t37 + t70 * t65;
t68 = (-t10 * t30 + t13 * t31) * r_i_i_C(1) + (-t10 * t31 - t13 * t30) * r_i_i_C(2);
t51 = t62 * t61;
t52 = t63 * t60;
t18 = t37 * t52 + (t37 * t51 + t65 * t58) * t33;
t23 = -t61 * t57 + t63 * t62;
t16 = t18 * t39 + t23 * t36;
t17 = -t65 * t52 + (t37 * t58 - t51 * t65) * t33;
t67 = (-t16 * t30 + t17 * t31) * r_i_i_C(1) + (-t16 * t31 - t17 * t30) * r_i_i_C(2);
t66 = r_i_i_C(3) + pkin(11) + pkin(10);
t38 = cos(qJ(5));
t48 = t38 * pkin(5) + r_i_i_C(1) * t31 - r_i_i_C(2) * t30 + pkin(4);
t35 = sin(qJ(5));
t46 = t35 * pkin(5) + t30 * r_i_i_C(1) + t31 * r_i_i_C(2) + pkin(9);
t45 = -t66 * t36 - t48 * t39 - pkin(3);
t1 = [0, t54, t45 * t13 + t46 * t14, t66 * t10 + t48 * (-t14 * t36 + t20 * t39) (-t10 * t35 + t13 * t38) * pkin(5) + t68, t68; 0, -t64, t45 * t11 + t46 * t12, t66 * t8 + t48 * (-t12 * t36 + t19 * t39) (t11 * t38 - t35 * t8) * pkin(5) + t69, t69; 1, t63, t45 * t17 + t46 * t18, t66 * t16 + t48 * (-t18 * t36 + t23 * t39) (-t16 * t35 + t17 * t38) * pkin(5) + t67, t67;];
Ja_transl  = t1;
