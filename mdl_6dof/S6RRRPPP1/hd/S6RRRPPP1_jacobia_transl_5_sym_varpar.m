% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:47
% EndTime: 2019-02-26 22:02:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (183->62), mult. (513->101), div. (0->0), fcn. (635->10), ass. (0->45)
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t24 = sin(pkin(6));
t47 = r_i_i_C(1) + qJ(4);
t65 = t47 * t24;
t67 = t30 * pkin(3) + t27 * t65 + pkin(2);
t28 = sin(qJ(2));
t26 = cos(pkin(6));
t41 = t47 * t26 + pkin(9);
t66 = t41 * t28;
t23 = sin(pkin(10));
t25 = cos(pkin(10));
t31 = cos(qJ(2));
t56 = t27 * t28;
t37 = t24 * t31 + t26 * t56;
t54 = t28 * t30;
t36 = t23 * t54 + t37 * t25;
t46 = r_i_i_C(3) + qJ(5);
t62 = r_i_i_C(2) - pkin(4);
t64 = t41 * t31 - t67 * t28 - t46 * t36 + t62 * (-t37 * t23 + t25 * t54);
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t49 = t32 * t27;
t51 = t30 * t31;
t19 = t29 * t51 - t49;
t48 = t32 * t30;
t52 = t29 * t27;
t18 = t31 * t52 + t48;
t55 = t28 * t29;
t40 = -t18 * t26 + t24 * t55;
t63 = -t19 * t23 + t40 * t25;
t20 = t29 * t30 - t31 * t49;
t60 = t20 * t24;
t59 = t23 * t26;
t58 = t25 * t26;
t57 = t26 * t30;
t53 = t28 * t32;
t50 = t31 * t26;
t45 = -t31 * pkin(2) - pkin(1);
t39 = t20 * t26 + t24 * t53;
t38 = -t24 * t28 + t27 * t50;
t21 = t31 * t48 + t52;
t15 = t26 * t53 - t60;
t3 = t21 * t23 - t39 * t25;
t1 = [-t19 * pkin(3) + t32 * pkin(8) - t62 * (-t19 * t25 - t40 * t23) - t18 * t65 + t46 * t63 + (t45 - t66) * t29, t64 * t32, t20 * pkin(3) - t62 * (t20 * t25 - t21 * t59) + t46 * (t20 * t23 + t21 * t58) + t21 * t65, t15, t3, 0; -qJ(4) * t60 + t21 * pkin(3) + t29 * pkin(8) + t15 * r_i_i_C(1) - t62 * (t21 * t25 + t39 * t23) + t46 * t3 + ((t26 * qJ(4) + pkin(9)) * t28 - t45) * t32, t64 * t29, -t18 * pkin(3) - t62 * (-t18 * t25 - t19 * t59) + t46 * (-t18 * t23 + t19 * t58) + t19 * t65, t18 * t24 + t26 * t55, -t63, 0; 0, -t62 * (-t38 * t23 + t25 * t51) + t46 * (t23 * t51 + t38 * t25) + t66 + t67 * t31 (t62 * (t23 * t57 + t25 * t27) - t46 * (t23 * t27 - t25 * t57) - pkin(3) * t27 + t30 * t65) * t28, t24 * t56 - t50, t36, 0;];
Ja_transl  = t1;
