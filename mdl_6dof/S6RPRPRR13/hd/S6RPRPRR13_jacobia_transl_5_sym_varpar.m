% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR13_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:38
% EndTime: 2019-02-26 20:55:38
% DurationCPUTime: 0.19s
% Computational Cost: add. (215->51), mult. (585->86), div. (0->0), fcn. (765->12), ass. (0->44)
t55 = pkin(4) + pkin(9);
t30 = cos(pkin(6));
t25 = sin(pkin(12));
t36 = cos(qJ(1));
t45 = t36 * t25;
t28 = cos(pkin(12));
t33 = sin(qJ(1));
t46 = t33 * t28;
t18 = t30 * t45 + t46;
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t43 = t36 * t28;
t48 = t33 * t25;
t17 = -t30 * t43 + t48;
t29 = cos(pkin(7));
t26 = sin(pkin(7));
t27 = sin(pkin(6));
t44 = t36 * t27;
t39 = t26 * t44;
t37 = t17 * t29 + t39;
t54 = -t18 * t35 + t32 * t37;
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t38 = t31 * r_i_i_C(1) + t34 * r_i_i_C(2) + qJ(4);
t53 = t18 * t32;
t51 = t26 * t30;
t50 = t27 * t28;
t49 = t29 * t35;
t47 = t33 * t27;
t42 = t27 * qJ(2);
t41 = r_i_i_C(3) + pkin(10) + pkin(3);
t40 = t26 * t47;
t11 = -t17 * t26 + t29 * t44;
t19 = -t30 * t46 - t45;
t13 = -t19 * t26 + t29 * t47;
t20 = -t30 * t48 + t43;
t16 = -t26 * t50 + t30 * t29;
t9 = t27 * t25 * t32 - t35 * t51 - t49 * t50;
t8 = t20 * t35 + (t19 * t29 + t40) * t32;
t7 = -t19 * t49 + t20 * t32 - t35 * t40;
t3 = t17 * t49 + t35 * t39 + t53;
t2 = t13 * t34 + t7 * t31;
t1 = -t13 * t31 + t7 * t34;
t4 = [-t18 * pkin(2) - t33 * pkin(1) + t36 * t42 + t41 * t54 + t38 * (-t35 * t37 - t53) + (t34 * r_i_i_C(1) - t31 * r_i_i_C(2) + t55) * t11, t47, t38 * t8 - t41 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t36 * pkin(1) + t20 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t55 * t13 + t33 * t42 + t41 * t8, -t44, -t41 * t3 - t38 * t54, t3 (t11 * t31 + t3 * t34) * r_i_i_C(1) + (t11 * t34 - t3 * t31) * r_i_i_C(2), 0; 0, t30, -t41 * t9 + t38 * (t32 * t51 + (t28 * t29 * t32 + t25 * t35) * t27) t9 (-t16 * t31 + t9 * t34) * r_i_i_C(1) + (-t16 * t34 - t9 * t31) * r_i_i_C(2), 0;];
Ja_transl  = t4;
