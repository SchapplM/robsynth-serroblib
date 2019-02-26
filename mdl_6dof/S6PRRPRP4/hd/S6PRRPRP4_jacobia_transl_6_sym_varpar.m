% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:04
% EndTime: 2019-02-26 20:03:04
% DurationCPUTime: 0.16s
% Computational Cost: add. (165->34), mult. (414->63), div. (0->0), fcn. (522->10), ass. (0->32)
t40 = pkin(5) + r_i_i_C(1);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t25 = t21 * r_i_i_C(2) + t40 * t18 + qJ(4);
t34 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(9);
t24 = -t25 * t19 - t34 * t22 - pkin(2);
t16 = sin(pkin(6));
t39 = t16 * t19;
t38 = t16 * t22;
t23 = cos(qJ(2));
t37 = t16 * t23;
t36 = cos(pkin(6));
t35 = cos(pkin(10));
t15 = sin(pkin(10));
t32 = t15 * t36;
t31 = t16 * t35;
t28 = t36 * t35;
t26 = -t18 * r_i_i_C(2) + t40 * t21 + pkin(4) + pkin(8);
t20 = sin(qJ(2));
t10 = t36 * t19 + t20 * t38;
t9 = t20 * t39 - t36 * t22;
t8 = -t20 * t32 + t35 * t23;
t7 = t35 * t20 + t23 * t32;
t6 = t15 * t23 + t20 * t28;
t5 = t15 * t20 - t23 * t28;
t4 = t15 * t39 + t8 * t22;
t3 = -t15 * t38 + t8 * t19;
t2 = -t19 * t31 + t6 * t22;
t1 = t6 * t19 + t22 * t31;
t11 = [0, t24 * t7 + t26 * t8, t25 * t4 - t34 * t3, t3 (-t3 * t18 - t7 * t21) * r_i_i_C(2) + t40 * (-t7 * t18 + t3 * t21) t4; 0, t24 * t5 + t26 * t6, -t34 * t1 + t25 * t2, t1 (-t1 * t18 - t5 * t21) * r_i_i_C(2) + t40 * (t1 * t21 - t5 * t18) t2; 1 (t26 * t20 - t24 * t23) * t16, t25 * t10 - t34 * t9, t9 (-t9 * t18 + t21 * t37) * r_i_i_C(2) + t40 * (t18 * t37 + t9 * t21) t10;];
Ja_transl  = t11;
