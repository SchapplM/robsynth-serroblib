% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:27
% EndTime: 2019-02-26 22:32:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (160->39), mult. (147->53), div. (0->0), fcn. (158->10), ass. (0->32)
t25 = cos(qJ(4));
t12 = pkin(4) * t25 + pkin(3);
t20 = qJ(2) + qJ(3);
t16 = sin(t20);
t17 = cos(t20);
t21 = -qJ(5) - pkin(9);
t43 = t17 * t12 + (r_i_i_C(3) - t21) * t16;
t18 = cos(qJ(2)) * pkin(2);
t42 = pkin(1) + t18 + t43;
t24 = sin(qJ(1));
t19 = qJ(4) + pkin(11);
t14 = sin(t19);
t38 = r_i_i_C(2) * t14;
t33 = t16 * t38;
t35 = t17 * t24;
t41 = r_i_i_C(3) * t35 + t24 * t33;
t22 = sin(qJ(4));
t40 = pkin(4) * t22;
t15 = cos(t19);
t39 = r_i_i_C(1) * t15;
t26 = cos(qJ(1));
t34 = t17 * t26;
t36 = r_i_i_C(3) * t34 + t26 * t33;
t32 = pkin(8) + pkin(7) + t40;
t30 = -t17 * t21 + (-t12 - t39) * t16;
t29 = (-t38 + t39) * t17 + t43;
t28 = -sin(qJ(2)) * pkin(2) + t30;
t4 = t14 * t24 + t15 * t34;
t3 = -t14 * t34 + t15 * t24;
t2 = t14 * t26 - t15 * t35;
t1 = t14 * t35 + t15 * t26;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t24 + t32 * t26, t28 * t26 + t36, t30 * t26 + t36, r_i_i_C(1) * t3 - r_i_i_C(2) * t4 + (-t22 * t34 + t24 * t25) * pkin(4), t26 * t16, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t32 * t24 + t42 * t26, t28 * t24 + t41, t30 * t24 + t41, -r_i_i_C(1) * t1 + t2 * r_i_i_C(2) + (-t22 * t35 - t25 * t26) * pkin(4), t24 * t16, 0; 0, t18 + t29, t29 (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15 - t40) * t16, -t17, 0;];
Ja_transl  = t5;
