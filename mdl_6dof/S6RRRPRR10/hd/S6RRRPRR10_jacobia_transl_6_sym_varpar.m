% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:23
% EndTime: 2019-02-26 22:21:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (227->48), mult. (405->73), div. (0->0), fcn. (490->10), ass. (0->36)
t28 = cos(qJ(2));
t24 = sin(qJ(2));
t39 = -r_i_i_C(3) - pkin(10) - pkin(9) + pkin(8);
t36 = t39 * t24;
t48 = t28 * pkin(2) + pkin(1) + t36;
t21 = qJ(5) + qJ(6);
t19 = sin(t21);
t20 = cos(t21);
t23 = sin(qJ(3));
t27 = cos(qJ(3));
t44 = (-(t19 * t27 - t20 * t23) * r_i_i_C(1) - (t19 * t23 + t20 * t27) * r_i_i_C(2)) * t24;
t22 = sin(qJ(5));
t37 = t22 * pkin(5) + qJ(4);
t34 = -t19 * r_i_i_C(1) - t37;
t32 = t20 * r_i_i_C(2) - t34;
t26 = cos(qJ(5));
t43 = t26 * pkin(5) + pkin(3) + pkin(4);
t35 = t19 * r_i_i_C(2) - t43;
t33 = t20 * r_i_i_C(1) - t35;
t47 = t32 * t23 + t33 * t27 + pkin(2);
t29 = cos(qJ(1));
t40 = t29 * t27;
t25 = sin(qJ(1));
t42 = t25 * t28;
t12 = t23 * t42 + t40;
t11 = t12 * t20;
t41 = t29 * t23;
t13 = t27 * t42 - t41;
t46 = (-t13 * t19 + t11) * r_i_i_C(1) + (-t12 * t19 - t13 * t20) * r_i_i_C(2);
t14 = -t25 * t27 + t28 * t41;
t15 = t25 * t23 + t28 * t40;
t5 = t14 * t20 - t15 * t19;
t6 = t14 * t19 + t15 * t20;
t45 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t31 = -t47 * t24 + t39 * t28;
t1 = [t29 * pkin(7) - t11 * r_i_i_C(2) + t34 * t12 - t33 * t13 - t48 * t25, t31 * t29, -t33 * t14 + t32 * t15, t14 (t14 * t26 - t15 * t22) * pkin(5) + t45, t45; t25 * pkin(7) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t37 * t14 + t43 * t15 + t48 * t29, t31 * t25, -t11 * r_i_i_C(1) + t35 * t12 + t32 * t13, t12 (t12 * t26 - t13 * t22) * pkin(5) + t46, t46; 0, t47 * t28 + t36 (-t43 * t23 + t37 * t27) * t24 - t44, t24 * t23 (-t22 * t27 + t23 * t26) * t24 * pkin(5) + t44, t44;];
Ja_transl  = t1;
