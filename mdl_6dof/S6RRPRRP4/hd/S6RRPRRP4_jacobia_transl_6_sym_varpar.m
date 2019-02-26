% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:45
% EndTime: 2019-02-26 21:47:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (241->35), mult. (210->45), div. (0->0), fcn. (241->10), ass. (0->34)
t22 = qJ(2) + pkin(10);
t17 = sin(t22);
t46 = r_i_i_C(2) + pkin(9) + pkin(8);
t51 = cos(qJ(2)) * pkin(2) + t46 * t17;
t48 = pkin(5) + r_i_i_C(1);
t39 = r_i_i_C(3) + qJ(6);
t28 = cos(qJ(4));
t15 = t28 * pkin(4) + pkin(3);
t18 = cos(t22);
t50 = t18 * t15 + pkin(1) + t51;
t23 = qJ(4) + qJ(5);
t19 = sin(t23);
t20 = cos(t23);
t49 = t39 * t19 + t48 * t20 + t15;
t25 = sin(qJ(4));
t47 = pkin(4) * t25;
t44 = t18 * t25;
t27 = sin(qJ(1));
t43 = t27 * t19;
t42 = t27 * t20;
t29 = cos(qJ(1));
t41 = t29 * t19;
t40 = t29 * t20;
t38 = t39 * t17 * t20;
t37 = t48 * t19;
t35 = qJ(3) + pkin(7) + t47;
t7 = t18 * t43 + t40;
t8 = t18 * t42 - t41;
t33 = t39 * t8 - t48 * t7;
t10 = t18 * t40 + t43;
t9 = t18 * t41 - t42;
t32 = t39 * t10 - t48 * t9;
t31 = -sin(qJ(2)) * pkin(2) + t46 * t18 - t49 * t17;
t1 = [-t50 * t27 + t35 * t29 - t39 * t7 - t48 * t8, t31 * t29, t27 (t27 * t28 - t29 * t44) * pkin(4) + t32, t32, t9; t48 * t10 + t35 * t27 + t50 * t29 + t39 * t9, t31 * t27, -t29 (-t27 * t44 - t28 * t29) * pkin(4) + t33, t33, t7; 0, t49 * t18 + t51, 0 (-t37 - t47) * t17 + t38, -t17 * t37 + t38, t17 * t19;];
Ja_transl  = t1;
