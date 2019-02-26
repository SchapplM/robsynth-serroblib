% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:37
% EndTime: 2019-02-26 20:00:37
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->36), mult. (380->65), div. (0->0), fcn. (481->12), ass. (0->32)
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t17 = pkin(11) + qJ(6);
t15 = sin(t17);
t16 = cos(t17);
t27 = pkin(5) * sin(pkin(11)) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + qJ(4);
t32 = pkin(3) + r_i_i_C(3) + pkin(9) + qJ(5);
t38 = t27 * t22 + t32 * t24 + pkin(2);
t20 = sin(pkin(6));
t37 = t20 * t22;
t36 = t20 * t24;
t25 = cos(qJ(2));
t35 = t20 * t25;
t34 = cos(pkin(6));
t33 = cos(pkin(10));
t19 = sin(pkin(10));
t31 = t19 * t34;
t30 = t20 * t33;
t29 = t34 * t33;
t28 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + pkin(8) + cos(pkin(11)) * pkin(5) + pkin(4);
t23 = sin(qJ(2));
t10 = t34 * t22 + t23 * t36;
t9 = t23 * t37 - t34 * t24;
t8 = -t23 * t31 + t33 * t25;
t7 = t33 * t23 + t25 * t31;
t6 = t19 * t25 + t23 * t29;
t5 = t19 * t23 - t25 * t29;
t4 = t19 * t37 + t8 * t24;
t3 = -t19 * t36 + t8 * t22;
t2 = -t22 * t30 + t6 * t24;
t1 = t6 * t22 + t24 * t30;
t11 = [0, t28 * t8 - t38 * t7, t27 * t4 - t32 * t3, t3, t4 (-t7 * t15 + t3 * t16) * r_i_i_C(1) + (-t3 * t15 - t7 * t16) * r_i_i_C(2); 0, t28 * t6 - t38 * t5, -t32 * t1 + t27 * t2, t1, t2 (t1 * t16 - t5 * t15) * r_i_i_C(1) + (-t1 * t15 - t5 * t16) * r_i_i_C(2); 1 (t28 * t23 + t38 * t25) * t20, t27 * t10 - t32 * t9, t9, t10 (t15 * t35 + t9 * t16) * r_i_i_C(1) + (-t9 * t15 + t16 * t35) * r_i_i_C(2);];
Ja_transl  = t11;
