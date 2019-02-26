% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (455->63), mult. (639->100), div. (0->0), fcn. (798->14), ass. (0->48)
t71 = r_i_i_C(3) + pkin(12) + pkin(11);
t49 = cos(qJ(5));
t36 = t49 * pkin(5) + pkin(4);
t42 = qJ(5) + qJ(6);
t38 = sin(t42);
t40 = cos(t42);
t60 = r_i_i_C(1) * t40 - r_i_i_C(2) * t38 + t36;
t50 = cos(qJ(3));
t37 = t50 * pkin(3) + pkin(2);
t43 = qJ(3) + qJ(4);
t39 = sin(t43);
t41 = cos(t43);
t76 = t71 * t39 + t60 * t41 + t37;
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t51 = cos(qJ(2));
t52 = cos(qJ(1));
t65 = cos(pkin(6));
t62 = t52 * t65;
t30 = t47 * t62 + t48 * t51;
t44 = sin(pkin(6));
t67 = t44 * t52;
t18 = t30 * t41 - t39 * t67;
t29 = t48 * t47 - t51 * t62;
t75 = (-t18 * t38 + t29 * t40) * r_i_i_C(1) + (-t18 * t40 - t29 * t38) * r_i_i_C(2);
t63 = t48 * t65;
t32 = -t47 * t63 + t52 * t51;
t69 = t44 * t48;
t22 = t32 * t41 + t39 * t69;
t31 = t52 * t47 + t51 * t63;
t5 = -t22 * t38 + t31 * t40;
t6 = t22 * t40 + t31 * t38;
t74 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t70 = t44 * t47;
t68 = t44 * t51;
t28 = t65 * t39 + t41 * t70;
t66 = (-t28 * t38 - t40 * t68) * r_i_i_C(1) + (-t28 * t40 + t38 * t68) * r_i_i_C(2);
t45 = sin(qJ(5));
t64 = t45 * pkin(5) + pkin(9) + pkin(10);
t46 = sin(qJ(3));
t61 = t44 * (pkin(3) * t46 + pkin(8));
t17 = -t30 * t39 - t41 * t67;
t59 = t60 * t17 + t71 * t18;
t21 = t32 * t39 - t41 * t69;
t58 = -t60 * t21 + t71 * t22;
t57 = t71 * t28 + t60 * (-t39 * t70 + t65 * t41);
t56 = t38 * r_i_i_C(1) + t40 * r_i_i_C(2) + t64;
t1 = [-t48 * pkin(1) + t71 * t17 - t60 * t18 - t56 * t29 - t30 * t37 + t52 * t61, -t31 * t76 + t56 * t32 (-t32 * t46 + t50 * t69) * pkin(3) + t58, t58 (-t22 * t45 + t31 * t49) * pkin(5) + t74, t74; t52 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t71 * t21 + t22 * t36 + t64 * t31 + t32 * t37 + t48 * t61, -t29 * t76 + t56 * t30 (-t30 * t46 - t50 * t67) * pkin(3) + t59, t59 (-t18 * t45 + t29 * t49) * pkin(5) + t75, t75; 0 (t56 * t47 + t76 * t51) * t44 (-t46 * t70 + t65 * t50) * pkin(3) + t57, t57 (-t28 * t45 - t49 * t68) * pkin(5) + t66, t66;];
Ja_transl  = t1;
