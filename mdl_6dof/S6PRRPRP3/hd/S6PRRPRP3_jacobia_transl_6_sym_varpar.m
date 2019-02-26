% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:22
% EndTime: 2019-02-26 20:02:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (262->47), mult. (487->84), div. (0->0), fcn. (623->12), ass. (0->41)
t28 = cos(pkin(11)) * pkin(4) + pkin(3);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t56 = r_i_i_C(2) + pkin(9) + qJ(4);
t58 = t28 * t38 + t56 * t36 + pkin(2);
t57 = pkin(5) + r_i_i_C(1);
t31 = pkin(11) + qJ(5);
t29 = sin(t31);
t55 = t29 * t38;
t30 = cos(t31);
t54 = t30 * t38;
t34 = sin(pkin(6));
t53 = t34 * t36;
t52 = t34 * t38;
t39 = cos(qJ(2));
t51 = t34 * t39;
t50 = t38 * t39;
t49 = r_i_i_C(3) + qJ(6);
t48 = cos(pkin(6));
t47 = cos(pkin(10));
t46 = pkin(4) * sin(pkin(11)) + pkin(8);
t33 = sin(pkin(10));
t44 = t33 * t48;
t43 = t34 * t47;
t42 = t48 * t47;
t40 = -t49 * t29 - t57 * t30 - t28;
t37 = sin(qJ(2));
t24 = t48 * t36 + t37 * t52;
t23 = t37 * t53 - t48 * t38;
t22 = -t37 * t44 + t47 * t39;
t21 = t47 * t37 + t39 * t44;
t20 = t33 * t39 + t37 * t42;
t19 = t33 * t37 - t39 * t42;
t14 = t22 * t38 + t33 * t53;
t13 = t22 * t36 - t33 * t52;
t12 = t20 * t38 - t36 * t43;
t11 = t20 * t36 + t38 * t43;
t9 = t24 * t29 + t30 * t51;
t3 = t14 * t29 - t21 * t30;
t1 = t12 * t29 - t19 * t30;
t2 = [0, t57 * (-t21 * t54 + t22 * t29) + t49 * (-t21 * t55 - t22 * t30) + t46 * t22 - t58 * t21, t40 * t13 + t56 * t14, t13, t49 * (t14 * t30 + t21 * t29) - t57 * t3, t3; 0, t57 * (-t19 * t54 + t20 * t29) + t49 * (-t19 * t55 - t20 * t30) + t46 * t20 - t58 * t19, t40 * t11 + t56 * t12, t11, t49 * (t12 * t30 + t19 * t29) - t57 * t1, t1; 1 (t57 * (t29 * t37 + t30 * t50) + t49 * (t29 * t50 - t30 * t37) + t46 * t37 + t58 * t39) * t34, t40 * t23 + t56 * t24, t23, -t57 * t9 + t49 * (t24 * t30 - t29 * t51) t9;];
Ja_transl  = t2;
