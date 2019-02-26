% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:49
% EndTime: 2019-02-26 22:47:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (306->43), mult. (209->56), div. (0->0), fcn. (221->12), ass. (0->43)
t29 = qJ(2) + qJ(3);
t25 = qJ(4) + t29;
t19 = sin(t25);
t20 = cos(t25);
t32 = cos(qJ(5));
t21 = t32 * pkin(5) + pkin(4);
t34 = -pkin(11) - pkin(10);
t60 = t20 * t21 + (r_i_i_C(3) - t34) * t19;
t28 = qJ(5) + qJ(6);
t24 = cos(t28);
t54 = r_i_i_C(1) * t24;
t39 = -t20 * t34 + (-t21 - t54) * t19;
t35 = t39 - pkin(3) * sin(t29);
t18 = pkin(3) * cos(t29);
t26 = cos(qJ(2)) * pkin(2);
t59 = pkin(1) + t18 + t26 + t60;
t33 = cos(qJ(1));
t45 = t33 * t24;
t22 = sin(t28);
t31 = sin(qJ(1));
t48 = t31 * t22;
t5 = t20 * t48 + t45;
t46 = t33 * t22;
t47 = t31 * t24;
t6 = -t20 * t47 + t46;
t58 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t20 * t46 + t47;
t8 = t20 * t45 + t48;
t57 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t30 = sin(qJ(5));
t55 = pkin(5) * t30;
t53 = r_i_i_C(2) * t22;
t43 = t19 * t53;
t50 = t20 * t31;
t51 = r_i_i_C(3) * t50 + t31 * t43;
t49 = t20 * t33;
t44 = r_i_i_C(3) * t49 + t33 * t43;
t42 = pkin(9) + pkin(8) + pkin(7) + t55;
t40 = -r_i_i_C(1) * t22 - r_i_i_C(2) * t24;
t38 = (-t53 + t54) * t20 + t60;
t37 = -sin(qJ(2)) * pkin(2) + t35;
t36 = t18 + t38;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t59 * t31 + t42 * t33, t37 * t33 + t44, t35 * t33 + t44, t39 * t33 + t44 (-t30 * t49 + t31 * t32) * pkin(5) + t57, t57; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t42 * t31 + t59 * t33, t37 * t31 + t51, t35 * t31 + t51, t39 * t31 + t51 (-t30 * t50 - t32 * t33) * pkin(5) + t58, t58; 0, t26 + t36, t36, t38 (t40 - t55) * t19, t40 * t19;];
Ja_transl  = t1;
