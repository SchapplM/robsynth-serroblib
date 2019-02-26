% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR15_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:47
% EndTime: 2019-02-26 22:38:47
% DurationCPUTime: 0.36s
% Computational Cost: add. (651->90), mult. (1797->159), div. (0->0), fcn. (2353->14), ass. (0->59)
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t54 = cos(qJ(1));
t75 = cos(pkin(6));
t67 = t54 * t75;
t81 = sin(qJ(1));
t37 = t81 * t50 - t53 * t67;
t38 = t50 * t67 + t81 * t53;
t49 = sin(qJ(3));
t74 = cos(pkin(7));
t68 = t49 * t74;
t45 = sin(pkin(7));
t46 = sin(pkin(6));
t77 = t46 * t54;
t71 = t45 * t77;
t82 = cos(qJ(3));
t16 = -t37 * t68 + t38 * t82 - t49 * t71;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t69 = t46 * t74;
t60 = t37 * t45 - t54 * t69;
t86 = t16 * t52 + t60 * t48;
t3 = t16 * t48 - t60 * t52;
t64 = t75 * t81;
t59 = t54 * t50 + t53 * t64;
t70 = t46 * t81;
t85 = -t45 * t70 + t59 * t74;
t84 = pkin(5) + pkin(11);
t83 = pkin(10) * t45;
t80 = t45 * t48;
t79 = t45 * t52;
t78 = t46 * t53;
t76 = t50 * t45;
t73 = r_i_i_C(3) + pkin(12) + pkin(4);
t72 = t46 * t76;
t66 = t75 * t45;
t63 = t74 * t82;
t47 = sin(qJ(6));
t51 = cos(qJ(6));
t62 = t47 * r_i_i_C(1) + t51 * r_i_i_C(2) + qJ(5);
t61 = t51 * r_i_i_C(1) - t47 * r_i_i_C(2) + t84;
t58 = -t45 * t78 + t75 * t74;
t56 = -t62 * t48 - t73 * t52 - pkin(3);
t15 = t37 * t63 + t38 * t49 + t82 * t71;
t55 = t59 * t45 + t81 * t69;
t39 = -t50 * t64 + t54 * t53;
t35 = (-t50 * t68 + t82 * t53) * t46;
t30 = t49 * t66 + (t82 * t50 + t53 * t68) * t46;
t29 = t46 * t50 * t49 - t63 * t78 - t82 * t66;
t24 = -t39 * t68 - t59 * t82;
t22 = -t37 * t82 - t38 * t68;
t20 = t39 * t82 - t85 * t49;
t19 = t39 * t49 + t85 * t82;
t13 = t30 * t48 - t52 * t58;
t8 = t20 * t52 + t48 * t55;
t7 = t20 * t48 - t52 * t55;
t2 = t19 * t51 + t7 * t47;
t1 = -t19 * t47 + t7 * t51;
t4 = [-t81 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) + pkin(9) * t77 - t60 * pkin(10) - t61 * t15 - t62 * t3 - t73 * t86, t24 * pkin(3) - t59 * pkin(2) + t39 * t83 + t62 * (t24 * t48 - t39 * t79) + t61 * (t39 * t63 - t59 * t49) + t73 * (t24 * t52 + t39 * t80) t56 * t19 + t61 * t20, t62 * t8 - t73 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t54 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + pkin(9) * t70 + t55 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t84 * t19 + t73 * t8, t38 * t83 - t37 * pkin(2) + t22 * pkin(3) + t62 * (t22 * t48 - t38 * t79) + t61 * (-t37 * t49 + t38 * t63) + t73 * (t22 * t52 + t38 * t80) t56 * t15 + t61 * t16, -t73 * t3 + t62 * t86, t3 (-t15 * t47 + t3 * t51) * r_i_i_C(1) + (-t15 * t51 - t3 * t47) * r_i_i_C(2); 0, t35 * pkin(3) + t62 * (t35 * t48 - t52 * t72) + t73 * (t35 * t52 + t48 * t72) + (t53 * pkin(2) + pkin(10) * t76 + t61 * (t49 * t53 + t50 * t63)) * t46, t56 * t29 + t61 * t30, t62 * (t30 * t52 + t48 * t58) - t73 * t13, t13 (t13 * t51 - t29 * t47) * r_i_i_C(1) + (-t13 * t47 - t29 * t51) * r_i_i_C(2);];
Ja_transl  = t4;
