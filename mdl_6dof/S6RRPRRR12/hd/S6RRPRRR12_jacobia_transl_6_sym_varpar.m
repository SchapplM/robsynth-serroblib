% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:23
% EndTime: 2019-02-26 22:00:24
% DurationCPUTime: 0.22s
% Computational Cost: add. (332->56), mult. (552->89), div. (0->0), fcn. (693->12), ass. (0->41)
t61 = r_i_i_C(3) + pkin(11);
t34 = sin(qJ(6));
t38 = cos(qJ(6));
t48 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
t58 = pkin(2) + pkin(10) + pkin(9);
t33 = sin(pkin(6));
t36 = sin(qJ(2));
t57 = t33 * t36;
t37 = sin(qJ(1));
t56 = t33 * t37;
t40 = cos(qJ(2));
t55 = t33 * t40;
t41 = cos(qJ(1));
t54 = t33 * t41;
t53 = cos(pkin(6));
t39 = cos(qJ(4));
t52 = t33 * (t39 * pkin(4) + pkin(3) + pkin(8));
t35 = sin(qJ(4));
t51 = t35 * pkin(4) + qJ(3);
t50 = t37 * t53;
t49 = t41 * t53;
t22 = t37 * t36 - t40 * t49;
t32 = qJ(4) + qJ(5);
t30 = sin(t32);
t31 = cos(t32);
t13 = -t22 * t30 + t31 * t54;
t11 = t22 * t31 + t30 * t54;
t47 = t34 * r_i_i_C(1) + t38 * r_i_i_C(2) + t58;
t24 = t41 * t36 + t40 * t50;
t10 = t24 * t30 + t31 * t56;
t9 = -t24 * t31 + t30 * t56;
t46 = t61 * t10 - t48 * t9;
t45 = t48 * t11 - t61 * t13;
t21 = -t30 * t55 + t53 * t31;
t44 = t61 * t21 + t48 * (-t53 * t30 - t31 * t55);
t43 = t48 * t30 - t61 * t31 + t51;
t25 = -t36 * t50 + t41 * t40;
t23 = t36 * t49 + t37 * t40;
t2 = t10 * t38 + t25 * t34;
t1 = -t10 * t34 + t25 * t38;
t3 = [-t37 * pkin(1) + t61 * t11 + t48 * t13 - t51 * t22 - t47 * t23 + t41 * t52, -t47 * t24 + t43 * t25, t24 (t24 * t39 - t35 * t56) * pkin(4) + t46, t46, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t41 * pkin(1) + t10 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t51 * t24 + t58 * t25 + t37 * t52 + t61 * t9, -t47 * t22 + t43 * t23, t22 (t22 * t39 + t35 * t54) * pkin(4) + t45, t45 (t13 * t34 + t23 * t38) * r_i_i_C(1) + (t13 * t38 - t23 * t34) * r_i_i_C(2); 0 (t43 * t36 + t47 * t40) * t33, -t55 (-t53 * t35 - t39 * t55) * pkin(4) + t44, t44 (-t21 * t34 + t38 * t57) * r_i_i_C(1) + (-t21 * t38 - t34 * t57) * r_i_i_C(2);];
Ja_transl  = t3;
