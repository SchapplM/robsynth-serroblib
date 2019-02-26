% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR10_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:49
% EndTime: 2019-02-26 21:19:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
t27 = cos(pkin(6));
t25 = cos(pkin(13));
t33 = cos(qJ(1));
t41 = t33 * t25;
t22 = sin(pkin(13));
t30 = sin(qJ(1));
t46 = t30 * t22;
t16 = -t27 * t41 + t46;
t23 = sin(pkin(7));
t26 = cos(pkin(7));
t24 = sin(pkin(6));
t42 = t33 * t24;
t11 = -t16 * t23 + t26 * t42;
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t43 = t33 * t22;
t44 = t30 * t25;
t17 = t27 * t43 + t44;
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t39 = t23 * t42;
t47 = t26 * t29;
t6 = t16 * t47 - t17 * t32 + t29 * t39;
t52 = t11 * t31 - t6 * t28;
t51 = t11 * t28 + t6 * t31;
t36 = t27 * t44 + t43;
t45 = t30 * t24;
t50 = -t23 * t45 + t36 * t26;
t49 = r_i_i_C(3) + pkin(10);
t48 = t23 * t27;
t40 = t24 * qJ(2);
t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
t34 = -t17 * t29 + (-t16 * t26 - t39) * t32;
t13 = t23 * t36 + t26 * t45;
t18 = -t27 * t46 + t41;
t15 = -t24 * t25 * t23 + t27 * t26;
t10 = t29 * t48 + (t22 * t32 + t25 * t47) * t24;
t8 = t18 * t32 - t50 * t29;
t7 = t18 * t29 + t50 * t32;
t2 = t13 * t28 + t8 * t31;
t1 = t13 * t31 - t8 * t28;
t3 = [-t30 * pkin(1) - t17 * pkin(2) + t6 * pkin(3) + t11 * pkin(9) + t51 * r_i_i_C(1) + t52 * r_i_i_C(2) + t33 * t40 + t49 * t34, t45, -t37 * t7 + t49 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t33 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t30 * t40 + t49 * t7, -t42, t37 * t34 - t49 * t6, -t52 * r_i_i_C(1) + t51 * r_i_i_C(2), 0, 0; 0, t27, t49 * t10 + t37 * (t32 * t48 + (t25 * t26 * t32 - t22 * t29) * t24) (-t10 * t28 + t15 * t31) * r_i_i_C(1) + (-t10 * t31 - t15 * t28) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
