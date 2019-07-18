% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobia_transl_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->28), mult. (107->38), div. (0->0), fcn. (117->8), ass. (0->30)
t16 = qJ(2) + qJ(4);
t14 = sin(t16);
t15 = cos(t16);
t18 = sin(qJ(5));
t38 = r_i_i_C(2) * t18;
t44 = pkin(4) + r_i_i_C(3);
t45 = t14 * t38 + t15 * t44;
t42 = t44 * t14;
t24 = pkin(2) + pkin(1);
t31 = cos(qJ(2)) * t24;
t41 = t31 + t42;
t21 = cos(qJ(5));
t39 = r_i_i_C(1) * t21;
t20 = sin(qJ(1));
t35 = t18 * t20;
t23 = cos(qJ(1));
t34 = t18 * t23;
t33 = t20 * t21;
t32 = t21 * t23;
t30 = t45 * t20;
t29 = t14 * t39;
t27 = t45 * t23;
t26 = -sin(qJ(2)) * t24 - t29;
t25 = (-t38 + t39) * t15 + t42;
t17 = -pkin(3) - qJ(3);
t4 = t15 * t32 + t35;
t3 = -t15 * t34 + t33;
t2 = -t15 * t33 + t34;
t1 = t15 * t35 + t32;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t17 - t41 * t20, t26 * t23 + t27, t20, -t23 * t29 + t27, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t17 + t41 * t23, t26 * t20 + t30, -t23, -t20 * t29 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t25 + t31, 0, t25, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t21) * t14;];
Ja_transl  = t5;
