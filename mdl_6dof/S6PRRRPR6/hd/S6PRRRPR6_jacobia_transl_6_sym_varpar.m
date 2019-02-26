% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (329->56), mult. (866->104), div. (0->0), fcn. (1129->12), ass. (0->40)
t31 = sin(qJ(3));
t35 = cos(qJ(3));
t45 = -r_i_i_C(3) - pkin(10) + pkin(9);
t54 = pkin(3) * t35 + t45 * t31 + pkin(2);
t28 = sin(pkin(6));
t53 = t28 * t31;
t52 = t28 * t35;
t36 = cos(qJ(2));
t51 = t28 * t36;
t30 = sin(qJ(4));
t50 = t30 * t35;
t34 = cos(qJ(4));
t49 = t34 * t35;
t48 = t35 * t36;
t47 = cos(pkin(6));
t46 = cos(pkin(11));
t27 = sin(pkin(11));
t43 = t27 * t47;
t42 = t28 * t46;
t41 = t47 * t46;
t29 = sin(qJ(6));
t33 = cos(qJ(6));
t40 = r_i_i_C(1) * t29 + r_i_i_C(2) * t33 + qJ(5);
t39 = r_i_i_C(1) * t33 - r_i_i_C(2) * t29 + pkin(4) + pkin(5);
t37 = t40 * t30 + t39 * t34 + pkin(3);
t32 = sin(qJ(2));
t24 = t47 * t31 + t32 * t52;
t22 = -t32 * t43 + t46 * t36;
t21 = t46 * t32 + t36 * t43;
t20 = t27 * t36 + t32 * t41;
t19 = t27 * t32 - t36 * t41;
t14 = t24 * t34 - t30 * t51;
t13 = t24 * t30 + t34 * t51;
t12 = t22 * t35 + t27 * t53;
t10 = t20 * t35 - t31 * t42;
t4 = t12 * t34 + t21 * t30;
t3 = t12 * t30 - t21 * t34;
t2 = t10 * t34 + t19 * t30;
t1 = t10 * t30 - t19 * t34;
t5 = [0, t22 * pkin(8) + t40 * (-t21 * t50 - t22 * t34) + t39 * (-t21 * t49 + t22 * t30) - t54 * t21, t45 * t12 + t37 * (-t22 * t31 + t27 * t52) -t39 * t3 + t40 * t4, t3 (-t29 * t4 + t3 * t33) * r_i_i_C(1) + (-t29 * t3 - t33 * t4) * r_i_i_C(2); 0, t20 * pkin(8) + t40 * (-t19 * t50 - t20 * t34) + t39 * (-t19 * t49 + t20 * t30) - t54 * t19, t45 * t10 + t37 * (-t20 * t31 - t35 * t42) -t39 * t1 + t40 * t2, t1 (t1 * t33 - t2 * t29) * r_i_i_C(1) + (-t1 * t29 - t2 * t33) * r_i_i_C(2); 1 (t40 * (t30 * t48 - t32 * t34) + t39 * (t30 * t32 + t34 * t48) + t32 * pkin(8) + t54 * t36) * t28, t45 * t24 + t37 * (-t32 * t53 + t47 * t35) -t39 * t13 + t40 * t14, t13 (t13 * t33 - t14 * t29) * r_i_i_C(1) + (-t13 * t29 - t14 * t33) * r_i_i_C(2);];
Ja_transl  = t5;
