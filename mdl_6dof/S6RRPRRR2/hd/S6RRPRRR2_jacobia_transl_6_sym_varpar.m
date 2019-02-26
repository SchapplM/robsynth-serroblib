% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:40
% EndTime: 2019-02-26 21:54:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (247->42), mult. (170->54), div. (0->0), fcn. (184->12), ass. (0->40)
t27 = qJ(2) + pkin(11);
t22 = qJ(4) + t27;
t19 = sin(t22);
t20 = cos(t22);
t31 = cos(qJ(5));
t21 = t31 * pkin(5) + pkin(4);
t33 = -pkin(10) - pkin(9);
t57 = t20 * t21 + (r_i_i_C(3) - t33) * t19;
t41 = pkin(3) * cos(t27) + cos(qJ(2)) * pkin(2);
t56 = pkin(1) + t41 + t57;
t28 = qJ(5) + qJ(6);
t24 = cos(t28);
t32 = cos(qJ(1));
t43 = t32 * t24;
t23 = sin(t28);
t30 = sin(qJ(1));
t46 = t30 * t23;
t5 = t20 * t46 + t43;
t44 = t32 * t23;
t45 = t30 * t24;
t6 = -t20 * t45 + t44;
t55 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t20 * t44 + t45;
t8 = t20 * t43 + t46;
t54 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t29 = sin(qJ(5));
t53 = pkin(5) * t29;
t52 = r_i_i_C(1) * t24;
t51 = r_i_i_C(2) * t23;
t40 = t19 * t51;
t48 = t20 * t30;
t49 = r_i_i_C(3) * t48 + t30 * t40;
t47 = t20 * t32;
t42 = r_i_i_C(3) * t47 + t32 * t40;
t39 = pkin(8) + qJ(3) + pkin(7) + t53;
t37 = -r_i_i_C(1) * t23 - r_i_i_C(2) * t24;
t36 = -t20 * t33 + (-t21 - t52) * t19;
t35 = (-t51 + t52) * t20 + t57;
t34 = -pkin(3) * sin(t27) - sin(qJ(2)) * pkin(2) + t36;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t56 * t30 + t39 * t32, t34 * t32 + t42, t30, t36 * t32 + t42 (-t29 * t47 + t30 * t31) * pkin(5) + t54, t54; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t39 * t30 + t56 * t32, t34 * t30 + t49, -t32, t36 * t30 + t49 (-t29 * t48 - t31 * t32) * pkin(5) + t55, t55; 0, t35 + t41, 0, t35 (t37 - t53) * t19, t37 * t19;];
Ja_transl  = t1;
