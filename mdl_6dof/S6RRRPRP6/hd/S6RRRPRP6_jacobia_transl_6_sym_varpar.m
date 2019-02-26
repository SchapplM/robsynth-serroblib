% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:03
% EndTime: 2019-02-26 22:12:03
% DurationCPUTime: 0.20s
% Computational Cost: add. (305->58), mult. (500->91), div. (0->0), fcn. (627->12), ass. (0->41)
t50 = pkin(5) + r_i_i_C(1);
t26 = -qJ(4) - pkin(9);
t27 = sin(qJ(5));
t31 = cos(qJ(5));
t36 = t31 * r_i_i_C(2) + t50 * t27 - t26;
t32 = cos(qJ(3));
t20 = t32 * pkin(3) + pkin(2);
t23 = qJ(3) + pkin(11);
t21 = sin(t23);
t22 = cos(t23);
t19 = t31 * pkin(5) + pkin(4);
t38 = t31 * r_i_i_C(1) - t27 * r_i_i_C(2) + t19;
t49 = r_i_i_C(3) + qJ(6) + pkin(10);
t51 = t49 * t21 + t38 * t22 + t20;
t24 = sin(pkin(6));
t29 = sin(qJ(2));
t48 = t24 * t29;
t30 = sin(qJ(1));
t47 = t24 * t30;
t33 = cos(qJ(2));
t46 = t24 * t33;
t34 = cos(qJ(1));
t45 = t24 * t34;
t44 = cos(pkin(6));
t41 = t34 * t44;
t12 = t29 * t41 + t30 * t33;
t4 = t12 * t22 - t21 * t45;
t42 = t30 * t44;
t28 = sin(qJ(3));
t40 = t24 * (pkin(3) * t28 + pkin(8));
t13 = t34 * t29 + t33 * t42;
t14 = -t29 * t42 + t34 * t33;
t8 = t14 * t22 + t21 * t47;
t1 = t13 * t31 - t8 * t27;
t3 = t12 * t21 + t22 * t45;
t11 = t30 * t29 - t33 * t41;
t10 = t44 * t21 + t22 * t48;
t9 = t21 * t48 - t44 * t22;
t7 = t14 * t21 - t22 * t47;
t2 = t13 * t27 + t8 * t31;
t5 = [-t30 * pkin(1) - t36 * t11 - t12 * t20 - t49 * t3 + t34 * t40 - t38 * t4, -t13 * t51 + t36 * t14, t49 * t8 + (-t14 * t28 + t32 * t47) * pkin(3) - t38 * t7, t13, -t2 * r_i_i_C(2) + t50 * t1, t7; t34 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t20 + t8 * t19 + t49 * t7 + t30 * t40 + (t27 * pkin(5) - t26) * t13, -t11 * t51 + t36 * t12, t49 * t4 + (-t12 * t28 - t32 * t45) * pkin(3) - t38 * t3, t11 (-t11 * t27 - t4 * t31) * r_i_i_C(2) + t50 * (t11 * t31 - t4 * t27) t3; 0 (t36 * t29 + t51 * t33) * t24, t49 * t10 + (-t28 * t48 + t44 * t32) * pkin(3) - t38 * t9, -t46 (-t10 * t31 + t27 * t46) * r_i_i_C(2) + t50 * (-t10 * t27 - t31 * t46) t9;];
Ja_transl  = t5;
