% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:57
% EndTime: 2019-02-26 21:56:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (278->43), mult. (377->63), div. (0->0), fcn. (449->10), ass. (0->35)
t26 = sin(qJ(2));
t30 = cos(qJ(2));
t48 = qJ(4) + qJ(5);
t45 = sin(t48);
t46 = cos(t48);
t16 = t26 * t45 + t30 * t46;
t27 = sin(qJ(1));
t10 = t16 * t27;
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t39 = t28 * r_i_i_C(1) - t24 * r_i_i_C(2) + pkin(5);
t53 = -r_i_i_C(3) - pkin(10);
t17 = t26 * t46 - t30 * t45;
t9 = t17 * t27;
t37 = -t53 * t10 + t39 * t9;
t31 = cos(qJ(1));
t40 = t31 * t45;
t41 = t31 * t46;
t11 = -t26 * t40 - t30 * t41;
t12 = -t26 * t41 + t30 * t40;
t36 = t53 * t11 - t39 * t12;
t35 = -t39 * t16 - t53 * t17;
t25 = sin(qJ(4));
t47 = pkin(4) * t25 + qJ(3);
t29 = cos(qJ(4));
t50 = t29 * pkin(4) + pkin(2) + pkin(3);
t34 = t47 * t26 + t50 * t30;
t54 = pkin(1) + t34;
t49 = pkin(7) - pkin(9) - pkin(8);
t44 = -t24 * r_i_i_C(1) - t28 * r_i_i_C(2);
t38 = pkin(4) * (-t25 * t30 + t26 * t29);
t33 = -t50 * t26 + t47 * t30;
t2 = -t11 * t28 - t27 * t24;
t1 = t11 * t24 - t27 * t28;
t3 = [-t53 * t9 - t39 * t10 + (t44 + t49) * t31 - t54 * t27, t33 * t31 - t36, t31 * t26, t31 * t38 + t36, t36, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t11 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t53 * t12 + t49 * t27 + t31 * t54, t33 * t27 - t37, t27 * t26, t27 * t38 + t37, t37 (-t10 * t24 + t31 * t28) * r_i_i_C(1) + (-t10 * t28 - t31 * t24) * r_i_i_C(2); 0, t34 - t35, -t30 (-t25 * t26 - t29 * t30) * pkin(4) + t35, t35, t44 * t17;];
Ja_transl  = t3;
