% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:37
% EndTime: 2019-02-26 20:06:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (154->35), mult. (303->65), div. (0->0), fcn. (382->12), ass. (0->32)
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t17 = pkin(12) + qJ(5);
t15 = sin(t17);
t16 = cos(t17);
t28 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + cos(pkin(12)) * pkin(4) + pkin(3);
t37 = r_i_i_C(3) + pkin(9) + qJ(4);
t38 = t37 * t22 + t28 * t24 + pkin(2);
t20 = sin(pkin(6));
t36 = t20 * t22;
t35 = t20 * t24;
t25 = cos(qJ(2));
t34 = t20 * t25;
t33 = cos(pkin(6));
t32 = cos(pkin(11));
t19 = sin(pkin(11));
t31 = t19 * t33;
t30 = t20 * t32;
t29 = t33 * t32;
t27 = sin(pkin(12)) * pkin(4) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
t23 = sin(qJ(2));
t10 = t33 * t22 + t23 * t35;
t9 = t23 * t36 - t33 * t24;
t8 = -t23 * t31 + t32 * t25;
t7 = t32 * t23 + t25 * t31;
t6 = t19 * t25 + t23 * t29;
t5 = t19 * t23 - t25 * t29;
t4 = t19 * t36 + t8 * t24;
t3 = -t19 * t35 + t8 * t22;
t2 = -t22 * t30 + t6 * t24;
t1 = t6 * t22 + t24 * t30;
t11 = [0, t27 * t8 - t38 * t7, -t28 * t3 + t37 * t4, t3 (-t4 * t15 + t7 * t16) * r_i_i_C(1) + (-t7 * t15 - t4 * t16) * r_i_i_C(2), 0; 0, t27 * t6 - t38 * t5, -t28 * t1 + t37 * t2, t1 (-t2 * t15 + t5 * t16) * r_i_i_C(1) + (-t5 * t15 - t2 * t16) * r_i_i_C(2), 0; 1 (t27 * t23 + t38 * t25) * t20, t37 * t10 - t28 * t9, t9 (-t10 * t15 - t16 * t34) * r_i_i_C(1) + (-t10 * t16 + t15 * t34) * r_i_i_C(2), 0;];
Ja_transl  = t11;
