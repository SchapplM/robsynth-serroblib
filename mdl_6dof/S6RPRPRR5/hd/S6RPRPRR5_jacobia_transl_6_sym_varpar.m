% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:15
% EndTime: 2019-02-26 20:51:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (226->38), mult. (273->54), div. (0->0), fcn. (332->9), ass. (0->31)
t17 = pkin(10) + qJ(3);
t15 = sin(t17);
t16 = cos(t17);
t36 = pkin(3) + pkin(4);
t23 = t15 * qJ(4) + t36 * t16;
t40 = cos(pkin(10)) * pkin(2) + pkin(1) + t23;
t32 = sin(qJ(5));
t33 = cos(qJ(5));
t8 = t15 * t33 - t16 * t32;
t19 = sin(qJ(6));
t21 = cos(qJ(6));
t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(5);
t20 = sin(qJ(1));
t3 = t8 * t20;
t35 = -r_i_i_C(3) - pkin(9);
t7 = t15 * t32 + t16 * t33;
t4 = t7 * t20;
t39 = t25 * t3 - t35 * t4;
t22 = cos(qJ(1));
t27 = t22 * t32;
t28 = t22 * t33;
t5 = -t15 * t27 - t16 * t28;
t6 = -t15 * t28 + t16 * t27;
t38 = -t25 * t6 + t35 * t5;
t37 = -t25 * t7 - t35 * t8;
t34 = -pkin(8) + pkin(7) + qJ(2);
t26 = -t19 * r_i_i_C(1) - t21 * r_i_i_C(2);
t24 = qJ(4) * t16 - t36 * t15;
t2 = -t20 * t19 - t5 * t21;
t1 = t5 * t19 - t20 * t21;
t9 = [-t35 * t3 - t25 * t4 + (t26 + t34) * t22 - t40 * t20, t20, t24 * t22 - t38, t22 * t15, t38, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t34 * t20 + t40 * t22 - t35 * t6, -t22, t24 * t20 - t39, t20 * t15, t39 (-t4 * t19 + t22 * t21) * r_i_i_C(1) + (-t22 * t19 - t4 * t21) * r_i_i_C(2); 0, 0, t23 - t37, -t16, t37, t26 * t8;];
Ja_transl  = t9;
