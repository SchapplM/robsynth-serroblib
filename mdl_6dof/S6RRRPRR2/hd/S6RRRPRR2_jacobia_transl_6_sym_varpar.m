% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:32
% EndTime: 2019-02-26 22:16:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (253->42), mult. (175->54), div. (0->0), fcn. (189->12), ass. (0->40)
t29 = qJ(2) + qJ(3);
t22 = pkin(11) + t29;
t18 = sin(t22);
t19 = cos(t22);
t32 = cos(qJ(5));
t21 = t32 * pkin(5) + pkin(4);
t34 = -pkin(10) - pkin(9);
t60 = (r_i_i_C(3) - t34) * t18 + t19 * t21 + pkin(3) * cos(t29);
t28 = qJ(5) + qJ(6);
t25 = cos(t28);
t53 = r_i_i_C(1) * t25;
t35 = (-t21 - t53) * t18 - t19 * t34 - pkin(3) * sin(t29);
t26 = cos(qJ(2)) * pkin(2);
t58 = pkin(1) + t26 + t60;
t33 = cos(qJ(1));
t44 = t33 * t25;
t23 = sin(t28);
t31 = sin(qJ(1));
t47 = t31 * t23;
t5 = t19 * t47 + t44;
t45 = t33 * t23;
t46 = t31 * t25;
t6 = -t19 * t46 + t45;
t57 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t19 * t45 + t46;
t8 = t19 * t44 + t47;
t56 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t30 = sin(qJ(5));
t54 = pkin(5) * t30;
t52 = r_i_i_C(2) * t23;
t42 = t18 * t52;
t49 = t19 * t31;
t50 = r_i_i_C(3) * t49 + t31 * t42;
t48 = t19 * t33;
t43 = r_i_i_C(3) * t48 + t33 * t42;
t41 = qJ(4) + pkin(8) + pkin(7) + t54;
t39 = -r_i_i_C(1) * t23 - r_i_i_C(2) * t25;
t37 = -sin(qJ(2)) * pkin(2) + t35;
t36 = (-t52 + t53) * t19 + t60;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t58 * t31 + t41 * t33, t37 * t33 + t43, t35 * t33 + t43, t31 (-t30 * t48 + t31 * t32) * pkin(5) + t56, t56; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t41 * t31 + t58 * t33, t37 * t31 + t50, t35 * t31 + t50, -t33 (-t30 * t49 - t32 * t33) * pkin(5) + t57, t57; 0, t26 + t36, t36, 0 (t39 - t54) * t18, t39 * t18;];
Ja_transl  = t1;
