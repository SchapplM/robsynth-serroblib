% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:09
% EndTime: 2019-02-26 22:01:09
% DurationCPUTime: 0.15s
% Computational Cost: add. (275->56), mult. (537->89), div. (0->0), fcn. (678->12), ass. (0->41)
t32 = sin(qJ(2));
t36 = cos(qJ(2));
t37 = cos(qJ(1));
t33 = sin(qJ(1));
t47 = cos(pkin(6));
t45 = t33 * t47;
t21 = t37 * t32 + t36 * t45;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t29 = sin(pkin(6));
t50 = t29 * t33;
t10 = t21 * t31 + t35 * t50;
t22 = -t32 * t45 + t37 * t36;
t28 = qJ(5) + qJ(6);
t26 = sin(t28);
t27 = cos(t28);
t5 = -t10 * t26 + t22 * t27;
t6 = t10 * t27 + t22 * t26;
t55 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t44 = t37 * t47;
t19 = t33 * t32 - t36 * t44;
t48 = t29 * t37;
t13 = -t19 * t31 + t35 * t48;
t20 = t32 * t44 + t33 * t36;
t54 = (t13 * t26 + t20 * t27) * r_i_i_C(1) + (t13 * t27 - t20 * t26) * r_i_i_C(2);
t49 = t29 * t36;
t18 = -t31 * t49 + t47 * t35;
t51 = t29 * t32;
t53 = (-t18 * t26 + t27 * t51) * r_i_i_C(1) + (-t18 * t27 - t26 * t51) * r_i_i_C(2);
t52 = r_i_i_C(3) + pkin(11) + pkin(10);
t46 = t29 * (pkin(3) + pkin(8));
t30 = sin(qJ(5));
t43 = t30 * pkin(5) + pkin(2) + pkin(9);
t34 = cos(qJ(5));
t25 = t34 * pkin(5) + pkin(4);
t42 = r_i_i_C(1) * t27 - r_i_i_C(2) * t26 + t25;
t41 = t19 * t35 + t31 * t48;
t40 = t26 * r_i_i_C(1) + t27 * r_i_i_C(2) + t43;
t39 = t42 * t31 - t52 * t35 + qJ(3);
t9 = -t21 * t35 + t31 * t50;
t1 = [-t33 * pkin(1) - t19 * qJ(3) + t42 * t13 - t40 * t20 + t37 * t46 + t52 * t41, -t40 * t21 + t39 * t22, t21, t52 * t10 - t42 * t9 (-t10 * t30 + t22 * t34) * pkin(5) + t55, t55; t37 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t21 * qJ(3) + t10 * t25 + t43 * t22 + t33 * t46 + t52 * t9, -t40 * t19 + t39 * t20, t19, -t52 * t13 + t42 * t41 (t13 * t30 + t20 * t34) * pkin(5) + t54, t54; 0 (t39 * t32 + t40 * t36) * t29, -t49, t52 * t18 + t42 * (-t47 * t31 - t35 * t49) (-t18 * t30 + t34 * t51) * pkin(5) + t53, t53;];
Ja_transl  = t1;
