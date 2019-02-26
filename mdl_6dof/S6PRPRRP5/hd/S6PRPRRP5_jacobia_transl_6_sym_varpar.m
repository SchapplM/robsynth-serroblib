% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:30
% EndTime: 2019-02-26 19:52:31
% DurationCPUTime: 0.17s
% Computational Cost: add. (145->36), mult. (354->63), div. (0->0), fcn. (447->10), ass. (0->32)
t40 = pkin(5) + r_i_i_C(1);
t20 = sin(qJ(5));
t23 = cos(qJ(5));
t41 = t23 * r_i_i_C(2) + t40 * t20 + pkin(2) + pkin(8);
t39 = r_i_i_C(3) + qJ(6) + pkin(9);
t16 = sin(pkin(6));
t21 = sin(qJ(4));
t38 = t16 * t21;
t22 = sin(qJ(2));
t37 = t16 * t22;
t24 = cos(qJ(4));
t36 = t16 * t24;
t25 = cos(qJ(2));
t35 = t16 * t25;
t18 = cos(pkin(6));
t34 = t18 * t22;
t33 = t18 * t25;
t29 = -t20 * r_i_i_C(2) + t40 * t23 + pkin(4);
t26 = t29 * t21 - t39 * t24 + qJ(3);
t17 = cos(pkin(10));
t15 = sin(pkin(10));
t12 = t18 * t24 - t21 * t35;
t11 = t18 * t21 + t24 * t35;
t10 = -t15 * t34 + t17 * t25;
t9 = t15 * t33 + t17 * t22;
t8 = t15 * t25 + t17 * t34;
t7 = t15 * t22 - t17 * t33;
t4 = t17 * t36 - t7 * t21;
t3 = t17 * t38 + t7 * t24;
t2 = t15 * t36 + t9 * t21;
t1 = t15 * t38 - t9 * t24;
t5 = [0, t26 * t10 - t41 * t9, t9, -t29 * t1 + t39 * t2 (-t10 * t20 - t2 * t23) * r_i_i_C(2) + t40 * (t10 * t23 - t2 * t20) t1; 0, t26 * t8 - t41 * t7, t7, t29 * t3 - t39 * t4 (-t8 * t20 + t4 * t23) * r_i_i_C(2) + t40 * (t4 * t20 + t8 * t23) -t3; 1 (t26 * t22 + t41 * t25) * t16, -t35, -t29 * t11 + t39 * t12 (-t12 * t23 - t20 * t37) * r_i_i_C(2) + t40 * (-t12 * t20 + t23 * t37) t11;];
Ja_transl  = t5;
