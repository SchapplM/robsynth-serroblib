% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR15_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:20
% EndTime: 2019-02-26 22:24:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (169->48), mult. (455->84), div. (0->0), fcn. (582->10), ass. (0->37)
t54 = pkin(10) + r_i_i_C(1);
t27 = sin(qJ(2));
t28 = sin(qJ(1));
t30 = cos(qJ(2));
t31 = cos(qJ(1));
t40 = cos(pkin(6));
t36 = t31 * t40;
t17 = t27 * t36 + t28 * t30;
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t16 = t28 * t27 - t30 * t36;
t23 = sin(pkin(7));
t25 = cos(pkin(7));
t24 = sin(pkin(6));
t48 = t24 * t31;
t33 = t16 * t25 + t23 * t48;
t53 = -t17 * t29 + t26 * t33;
t1 = t17 * t26 + t29 * t33;
t52 = -r_i_i_C(2) + pkin(3);
t49 = t24 * t28;
t47 = t25 * t26;
t46 = t25 * t29;
t45 = t26 * t27;
t44 = t26 * t30;
t43 = t27 * t29;
t42 = t29 * t30;
t41 = r_i_i_C(3) + qJ(4);
t39 = t23 * t49;
t38 = t54 * t23;
t37 = t28 * t40;
t35 = t40 * t23;
t19 = -t27 * t37 + t31 * t30;
t18 = -t31 * t27 - t30 * t37;
t11 = -t29 * t35 + (-t25 * t42 + t45) * t24;
t6 = t19 * t29 + (t18 * t25 + t39) * t26;
t5 = -t18 * t46 + t19 * t26 - t29 * t39;
t2 = [-t17 * pkin(2) - t28 * pkin(1) + pkin(9) * t48 + t52 * t53 - t41 * t1 + t54 * (-t16 * t23 + t25 * t48) t18 * pkin(2) + t41 * (t18 * t26 + t19 * t46) + t19 * t38 + t52 * (t18 * t29 - t19 * t47) t41 * t6 - t5 * t52, t5, 0, 0; t31 * pkin(1) + t19 * pkin(2) + pkin(9) * t49 + t41 * t5 + t52 * t6 + t54 * (-t18 * t23 + t25 * t49) -t16 * pkin(2) + t52 * (-t16 * t29 - t17 * t47) + t41 * (-t16 * t26 + t17 * t46) + t17 * t38, -t52 * t1 - t41 * t53, t1, 0, 0; 0 (t52 * (-t25 * t45 + t42) + t41 * (t25 * t43 + t44) + pkin(2) * t30 + t27 * t38) * t24, t41 * (t26 * t35 + (t25 * t44 + t43) * t24) - t52 * t11, t11, 0, 0;];
Ja_transl  = t2;
