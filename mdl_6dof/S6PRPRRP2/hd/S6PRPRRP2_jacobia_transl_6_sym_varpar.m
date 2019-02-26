% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:04
% EndTime: 2019-02-26 19:51:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (330->52), mult. (857->93), div. (0->0), fcn. (1136->12), ass. (0->39)
t37 = sin(pkin(11));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t53 = cos(pkin(11));
t50 = -t44 * t37 + t47 * t53;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t61 = pkin(9) + r_i_i_C(2);
t49 = pkin(4) * t46 + t61 * t43 + pkin(3);
t62 = pkin(5) + r_i_i_C(1);
t38 = sin(pkin(10));
t39 = sin(pkin(6));
t60 = t38 * t39;
t40 = cos(pkin(10));
t59 = t40 * t39;
t41 = cos(pkin(6));
t58 = t41 * t47;
t42 = sin(qJ(5));
t57 = t42 * t46;
t45 = cos(qJ(5));
t55 = t45 * t46;
t54 = r_i_i_C(3) + qJ(6);
t34 = -t47 * t37 - t44 * t53;
t32 = t34 * t41;
t21 = -t40 * t32 + t38 * t50;
t22 = -t38 * t32 - t40 * t50;
t48 = t54 * t42 + t62 * t45 + pkin(4);
t31 = t50 * t41;
t30 = t34 * t39;
t29 = t50 * t39;
t26 = -t30 * t46 + t41 * t43;
t23 = -t38 * t31 + t40 * t34;
t20 = t40 * t31 + t38 * t34;
t14 = -t22 * t46 + t43 * t60;
t12 = t21 * t46 - t43 * t59;
t9 = t26 * t42 + t29 * t45;
t3 = t14 * t42 + t23 * t45;
t1 = t12 * t42 + t20 * t45;
t2 = [0, -t22 * pkin(8) + t62 * (-t22 * t42 + t23 * t55) + t54 * (t22 * t45 + t23 * t57) + (-t38 * t58 - t40 * t44) * pkin(2) + t49 * t23, t60, t61 * t14 + t48 * (t22 * t43 + t46 * t60) t54 * (t14 * t45 - t23 * t42) - t62 * t3, t3; 0, t21 * pkin(8) + t62 * (t20 * t55 + t21 * t42) + t54 * (t20 * t57 - t21 * t45) + (-t38 * t44 + t40 * t58) * pkin(2) + t49 * t20, -t59, t61 * t12 + t48 * (-t21 * t43 - t46 * t59) t54 * (t12 * t45 - t20 * t42) - t62 * t1, t1; 1, t39 * t47 * pkin(2) - t30 * pkin(8) + t62 * (t29 * t55 - t30 * t42) + t54 * (t29 * t57 + t30 * t45) + t49 * t29, t41, t61 * t26 + t48 * (t30 * t43 + t41 * t46) -t62 * t9 + t54 * (t26 * t45 - t29 * t42) t9;];
Ja_transl  = t2;
