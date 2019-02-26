% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (238->34), mult. (205->44), div. (0->0), fcn. (236->9), ass. (0->34)
t46 = pkin(5) + r_i_i_C(1);
t37 = r_i_i_C(3) + qJ(6);
t26 = cos(qJ(4));
t16 = t26 * pkin(4) + pkin(3);
t21 = pkin(10) + qJ(3);
t18 = cos(t21);
t17 = sin(t21);
t44 = r_i_i_C(2) + pkin(9) + pkin(8);
t34 = t44 * t17;
t48 = t34 + t18 * t16 + cos(pkin(10)) * pkin(2) + pkin(1);
t22 = qJ(4) + qJ(5);
t19 = sin(t22);
t20 = cos(t22);
t47 = t37 * t19 + t46 * t20 + t16;
t24 = sin(qJ(4));
t45 = pkin(4) * t24;
t42 = t18 * t24;
t25 = sin(qJ(1));
t41 = t25 * t19;
t40 = t25 * t20;
t27 = cos(qJ(1));
t39 = t27 * t19;
t38 = t27 * t20;
t36 = t37 * t17 * t20;
t35 = t46 * t19;
t33 = pkin(7) + qJ(2) + t45;
t7 = t18 * t41 + t38;
t8 = t18 * t40 - t39;
t31 = t37 * t8 - t46 * t7;
t10 = t18 * t38 + t41;
t9 = t18 * t39 - t40;
t30 = t37 * t10 - t46 * t9;
t29 = -t47 * t17 + t44 * t18;
t1 = [-t48 * t25 + t33 * t27 - t37 * t7 - t46 * t8, t25, t29 * t27 (t25 * t26 - t27 * t42) * pkin(4) + t30, t30, t9; t46 * t10 + t33 * t25 + t48 * t27 + t37 * t9, -t27, t29 * t25 (-t25 * t42 - t26 * t27) * pkin(4) + t31, t31, t7; 0, 0, t47 * t18 + t34 (-t35 - t45) * t17 + t36, -t17 * t35 + t36, t17 * t19;];
Ja_transl  = t1;
