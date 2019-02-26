% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:42
% EndTime: 2019-02-26 22:13:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (194->43), mult. (477->62), div. (0->0), fcn. (596->8), ass. (0->35)
t31 = sin(qJ(5));
t32 = sin(qJ(3));
t35 = cos(qJ(5));
t36 = cos(qJ(3));
t51 = r_i_i_C(3) + qJ(6);
t54 = t35 * t36;
t57 = t31 * t36;
t59 = pkin(5) + r_i_i_C(1);
t60 = pkin(3) + pkin(4);
t71 = -qJ(4) * t32 - t60 * t36 - pkin(2) - t59 * (t31 * t32 + t54) - t51 * (-t32 * t35 + t57);
t37 = cos(qJ(2));
t33 = sin(qJ(2));
t50 = pkin(8) - pkin(9) - r_i_i_C(2);
t45 = t50 * t33;
t70 = t37 * pkin(2) + pkin(1) + t45;
t38 = cos(qJ(1));
t52 = t38 * t36;
t34 = sin(qJ(1));
t55 = t34 * t37;
t22 = t32 * t55 + t52;
t53 = t38 * t32;
t23 = t36 * t55 - t53;
t44 = t22 * t31 + t23 * t35;
t65 = -t22 * t35 + t23 * t31;
t69 = t51 * t44 - t59 * t65;
t64 = t71 * t33 + t50 * t37;
t24 = -t34 * t36 + t37 * t53;
t25 = t34 * t32 + t37 * t52;
t10 = t24 * t31 + t25 * t35;
t9 = -t24 * t35 + t25 * t31;
t63 = t51 * t10 - t59 * t9;
t56 = t33 * t32;
t19 = t33 * t57 - t35 * t56;
t61 = -t51 * (-t31 * t56 - t33 * t54) - t59 * t19;
t1 = [t38 * pkin(7) - t22 * qJ(4) - t60 * t23 - t70 * t34 - t59 * t44 - t51 * t65, t64 * t38, t25 * qJ(4) - t60 * t24 - t63, t24, t63, t9; t34 * pkin(7) + t24 * qJ(4) + t59 * t10 + t60 * t25 + t70 * t38 + t51 * t9, t64 * t34, t23 * qJ(4) - t60 * t22 - t69, t22, t69, t65; 0, -t71 * t37 + t45 (qJ(4) * t36 - t60 * t32) * t33 - t61, t56, t61, t19;];
Ja_transl  = t1;
