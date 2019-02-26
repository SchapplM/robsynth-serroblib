% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP3
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
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:02
% EndTime: 2019-02-26 21:09:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (252->34), mult. (205->47), div. (0->0), fcn. (234->10), ass. (0->32)
t42 = pkin(5) + r_i_i_C(1);
t35 = r_i_i_C(3) + qJ(6);
t24 = cos(qJ(4));
t15 = t24 * pkin(4) + pkin(3);
t25 = cos(qJ(3));
t23 = sin(qJ(3));
t40 = r_i_i_C(2) + pkin(9) + pkin(8);
t31 = t40 * t23;
t44 = t25 * t15 + pkin(2) + t31;
t21 = qJ(4) + qJ(5);
t18 = sin(t21);
t19 = cos(t21);
t43 = t35 * t18 + t42 * t19 + t15;
t22 = sin(qJ(4));
t41 = pkin(4) * t22;
t39 = t18 * t25;
t37 = t19 * t25;
t36 = t22 * t25;
t34 = t35 * t19 * t23;
t33 = t42 * t18;
t32 = pkin(7) + t41;
t20 = qJ(1) + pkin(10);
t16 = sin(t20);
t17 = cos(t20);
t7 = t16 * t39 + t17 * t19;
t8 = t16 * t37 - t17 * t18;
t29 = t35 * t8 - t42 * t7;
t10 = t16 * t18 + t17 * t37;
t9 = -t16 * t19 + t17 * t39;
t28 = t35 * t10 - t42 * t9;
t27 = -t43 * t23 + t40 * t25;
t1 = [-sin(qJ(1)) * pkin(1) - t42 * t8 - t35 * t7 + t32 * t17 - t44 * t16, 0, t27 * t17 (t16 * t24 - t17 * t36) * pkin(4) + t28, t28, t9; cos(qJ(1)) * pkin(1) + t35 * t9 + t32 * t16 + t42 * t10 + t44 * t17, 0, t27 * t16 (-t16 * t36 - t17 * t24) * pkin(4) + t29, t29, t7; 0, 1, t43 * t25 + t31 (-t33 - t41) * t23 + t34, -t23 * t33 + t34, t23 * t18;];
Ja_transl  = t1;
