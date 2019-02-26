% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR14_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:11
% EndTime: 2019-02-26 22:38:11
% DurationCPUTime: 0.38s
% Computational Cost: add. (488->100), mult. (1351->170), div. (0->0), fcn. (1766->14), ass. (0->64)
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t52 = cos(qJ(1));
t67 = cos(pkin(6));
t63 = t52 * t67;
t81 = sin(qJ(1));
t34 = t81 * t48 - t51 * t63;
t35 = t48 * t63 + t81 * t51;
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t42 = sin(pkin(7));
t43 = sin(pkin(6));
t76 = t43 * t52;
t65 = t42 * t76;
t45 = cos(pkin(7));
t75 = t45 * t47;
t16 = t34 * t75 - t35 * t50 + t47 * t65;
t29 = -t34 * t42 + t45 * t76;
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t83 = t16 * t46 - t29 * t49;
t4 = t16 * t49 + t29 * t46;
t41 = sin(pkin(13));
t44 = cos(pkin(13));
t59 = r_i_i_C(1) * t44 - r_i_i_C(2) * t41 + pkin(4);
t68 = r_i_i_C(3) + qJ(5);
t54 = t68 * t46 + t59 * t49 + pkin(3);
t82 = pkin(10) * t42;
t78 = t42 * t46;
t77 = t42 * t49;
t74 = t45 * t50;
t73 = t47 * t48;
t72 = t47 * t51;
t71 = t48 * t42;
t70 = t48 * t50;
t69 = t50 * t51;
t66 = t43 * t71;
t64 = t43 * t81;
t62 = t67 * t42;
t61 = t42 * t64;
t60 = t67 * t81;
t58 = t41 * r_i_i_C(1) + t44 * r_i_i_C(2) + pkin(11);
t57 = -t43 * t51 * t42 + t67 * t45;
t56 = t52 * t48 + t51 * t60;
t55 = t56 * t50;
t15 = -t35 * t47 + (-t34 * t45 - t65) * t50;
t53 = t56 * t42 + t45 * t64;
t36 = -t48 * t60 + t52 * t51;
t32 = (-t45 * t73 + t69) * t43;
t31 = (t45 * t70 + t72) * t43;
t28 = t47 * t62 + (t45 * t72 + t70) * t43;
t24 = t32 * t49 + t46 * t66;
t22 = -t36 * t75 - t55;
t21 = t36 * t74 - t56 * t47;
t20 = -t34 * t50 - t35 * t75;
t19 = -t34 * t47 + t35 * t74;
t18 = t36 * t50 + (-t56 * t45 + t61) * t47;
t17 = t36 * t47 + t45 * t55 - t50 * t61;
t11 = t28 * t46 - t57 * t49;
t10 = t22 * t49 + t36 * t78;
t8 = t20 * t49 + t35 * t78;
t6 = t18 * t49 + t53 * t46;
t5 = t18 * t46 - t53 * t49;
t1 = [(t15 * t41 + t4 * t44) * r_i_i_C(1) + (t15 * t44 - t4 * t41) * r_i_i_C(2) + t4 * pkin(4) + t16 * pkin(3) + t15 * pkin(11) - t35 * pkin(2) - t81 * pkin(1) + pkin(9) * t76 + t68 * t83 + t29 * pkin(10) (t10 * t44 + t21 * t41) * r_i_i_C(1) + (-t10 * t41 + t21 * t44) * r_i_i_C(2) + t10 * pkin(4) + t22 * pkin(3) + t21 * pkin(11) - t56 * pkin(2) + t36 * t82 + t68 * (t22 * t46 - t36 * t77) -t54 * t17 + t58 * t18, -t59 * t5 + t68 * t6, t5, 0; (t17 * t41 + t6 * t44) * r_i_i_C(1) + (t17 * t44 - t6 * t41) * r_i_i_C(2) + t6 * pkin(4) + t18 * pkin(3) + t17 * pkin(11) + t36 * pkin(2) + t52 * pkin(1) + pkin(9) * t64 + t68 * t5 + t53 * pkin(10) (t19 * t41 + t8 * t44) * r_i_i_C(1) + (t19 * t44 - t8 * t41) * r_i_i_C(2) + t8 * pkin(4) + t20 * pkin(3) + t19 * pkin(11) - t34 * pkin(2) + t35 * t82 + t68 * (t20 * t46 - t35 * t77) t54 * t15 - t16 * t58, -t68 * t4 + t59 * t83, -t83, 0; 0 (t24 * t44 + t31 * t41) * r_i_i_C(1) + (-t24 * t41 + t31 * t44) * r_i_i_C(2) + t24 * pkin(4) + t32 * pkin(3) + t31 * pkin(11) + (t51 * pkin(2) + pkin(10) * t71) * t43 + t68 * (t32 * t46 - t49 * t66) t58 * t28 + t54 * (t50 * t62 + (t45 * t69 - t73) * t43) t68 * (t28 * t49 + t57 * t46) - t59 * t11, t11, 0;];
Ja_transl  = t1;
