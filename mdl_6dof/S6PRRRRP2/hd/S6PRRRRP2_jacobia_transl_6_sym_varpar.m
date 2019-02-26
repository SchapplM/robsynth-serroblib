% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:42
% EndTime: 2019-02-26 20:15:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (376->53), mult. (615->90), div. (0->0), fcn. (778->12), ass. (0->44)
t69 = r_i_i_C(3) + qJ(6);
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t83 = pkin(5) + r_i_i_C(1);
t88 = t69 * t52 + t55 * t83 + pkin(4);
t87 = pkin(10) + r_i_i_C(2);
t48 = qJ(3) + qJ(4);
t46 = sin(t48);
t47 = cos(t48);
t56 = cos(qJ(3));
t84 = t56 * pkin(3) + pkin(4) * t47 + t87 * t46 + pkin(2);
t77 = t47 * t52;
t76 = t47 * t55;
t49 = sin(pkin(11));
t50 = sin(pkin(6));
t75 = t49 * t50;
t54 = sin(qJ(2));
t74 = t49 * t54;
t57 = cos(qJ(2));
t73 = t49 * t57;
t72 = t50 * t54;
t71 = t52 * t57;
t70 = t55 * t57;
t67 = cos(pkin(11));
t65 = t50 * t67;
t64 = t67 * t54;
t63 = t67 * t57;
t51 = cos(pkin(6));
t40 = t51 * t64 + t73;
t22 = t40 * t47 - t46 * t65;
t61 = t87 * t22 + t88 * (-t40 * t46 - t47 * t65);
t42 = -t51 * t74 + t63;
t24 = t42 * t47 + t46 * t75;
t60 = t87 * t24 + t88 * (-t42 * t46 + t47 * t75);
t36 = t51 * t46 + t47 * t72;
t59 = t87 * t36 + t88 * (-t46 * t72 + t51 * t47);
t58 = -pkin(9) - pkin(8);
t53 = sin(qJ(3));
t41 = t51 * t73 + t64;
t39 = -t51 * t63 + t74;
t25 = t36 * t52 + t50 * t70;
t3 = t24 * t52 - t41 * t55;
t1 = t22 * t52 - t39 * t55;
t2 = [0, -t42 * t58 + t83 * (-t41 * t76 + t42 * t52) + t69 * (-t41 * t77 - t42 * t55) - t84 * t41 (-t42 * t53 + t56 * t75) * pkin(3) + t60, t60, t69 * (t24 * t55 + t41 * t52) - t83 * t3, t3; 0, -t40 * t58 + t83 * (-t39 * t76 + t40 * t52) + t69 * (-t39 * t77 - t40 * t55) - t84 * t39 (-t40 * t53 - t56 * t65) * pkin(3) + t61, t61, t69 * (t22 * t55 + t39 * t52) - t83 * t1, t1; 1 (t83 * (t47 * t70 + t52 * t54) + t69 * (t47 * t71 - t54 * t55) - t54 * t58 + t84 * t57) * t50 (t51 * t56 - t53 * t72) * pkin(3) + t59, t59, t69 * (t36 * t55 - t50 * t71) - t83 * t25, t25;];
Ja_transl  = t2;
