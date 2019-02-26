% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:53
% EndTime: 2019-02-26 22:24:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (252->39), mult. (222->51), div. (0->0), fcn. (246->10), ass. (0->36)
t22 = qJ(4) + pkin(10);
t17 = sin(t22);
t18 = cos(t22);
t36 = r_i_i_C(3) + qJ(6);
t47 = pkin(5) + r_i_i_C(1);
t51 = t36 * t17 + t47 * t18;
t28 = cos(qJ(4));
t15 = t28 * pkin(4) + pkin(3);
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t20 = cos(t23);
t24 = -qJ(5) - pkin(9);
t49 = t20 * t15 + (r_i_i_C(2) - t24) * t19;
t21 = cos(qJ(2)) * pkin(2);
t48 = pkin(1) + t21 + t49;
t25 = sin(qJ(4));
t46 = pkin(4) * t25;
t27 = sin(qJ(1));
t42 = t20 * t27;
t29 = cos(qJ(1));
t41 = t20 * t29;
t40 = t27 * t17;
t39 = t27 * t18;
t38 = t29 * t17;
t37 = t29 * t18;
t35 = pkin(8) + pkin(7) + t46;
t33 = t51 * t20 + t49;
t32 = -t20 * t24 + (-t15 - t51) * t19;
t31 = -sin(qJ(2)) * pkin(2) + t32;
t11 = r_i_i_C(2) * t41;
t10 = r_i_i_C(2) * t42;
t4 = t20 * t37 + t40;
t3 = t20 * t38 - t39;
t2 = t20 * t39 - t38;
t1 = t20 * t40 + t37;
t5 = [-t36 * t1 - t47 * t2 - t48 * t27 + t35 * t29, t31 * t29 + t11, t32 * t29 + t11, t36 * t4 - t47 * t3 + (-t25 * t41 + t27 * t28) * pkin(4), t29 * t19, t3; t35 * t27 + t48 * t29 + t36 * t3 + t47 * t4, t31 * t27 + t10, t32 * t27 + t10, t36 * t2 - t47 * t1 + (-t25 * t42 - t28 * t29) * pkin(4), t27 * t19, t1; 0, t21 + t33, t33 (-t47 * t17 + t36 * t18 - t46) * t19, -t20, t19 * t17;];
Ja_transl  = t5;
