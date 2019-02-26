% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP4
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

function Ja_transl = S6RPRRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:48
% EndTime: 2019-02-26 21:09:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (192->39), mult. (149->45), div. (0->0), fcn. (162->9), ass. (0->32)
t20 = pkin(10) + qJ(3);
t18 = qJ(4) + t20;
t14 = sin(t18);
t15 = cos(t18);
t22 = sin(qJ(5));
t39 = r_i_i_C(2) * t22;
t45 = r_i_i_C(3) * t15 + t14 * t39;
t24 = cos(qJ(5));
t16 = pkin(5) * t24 + pkin(4);
t21 = -qJ(6) - pkin(9);
t44 = t15 * t16 + (r_i_i_C(3) - t21) * t14;
t43 = pkin(5) + r_i_i_C(1);
t13 = pkin(3) * cos(t20);
t42 = t13 + cos(pkin(10)) * pkin(2) + pkin(1) + t44;
t23 = sin(qJ(1));
t41 = t45 * t23;
t40 = r_i_i_C(1) * t24;
t25 = cos(qJ(1));
t36 = t45 * t25;
t35 = t22 * t25;
t34 = t23 * t22;
t33 = t23 * t24;
t32 = t24 * t25;
t30 = pkin(5) * t22 + pkin(7) + pkin(8) + qJ(2);
t3 = -t15 * t35 + t33;
t1 = t15 * t34 + t32;
t28 = -t15 * t21 + (-t16 - t40) * t14;
t27 = (-t39 + t40) * t15 + t44;
t26 = -pkin(3) * sin(t20) + t28;
t4 = t15 * t32 + t34;
t2 = -t15 * t33 + t35;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t23 + t30 * t25, t23, t26 * t25 + t36, t28 * t25 + t36, -t4 * r_i_i_C(2) + t43 * t3, t25 * t14; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t30 * t23 + t42 * t25, -t25, t26 * t23 + t41, t28 * t23 + t41, t2 * r_i_i_C(2) - t43 * t1, t23 * t14; 0, 0, t13 + t27, t27 (-r_i_i_C(2) * t24 - t43 * t22) * t14, -t15;];
Ja_transl  = t5;
