% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP5
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

function Ja_transl = S6PRRPRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:30
% EndTime: 2019-02-26 20:03:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (125->33), mult. (325->63), div. (0->0), fcn. (411->10), ass. (0->28)
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t24 = t16 * r_i_i_C(1) + t19 * r_i_i_C(2) + qJ(4);
t28 = pkin(3) + pkin(9) + r_i_i_C(3);
t34 = t24 * t17 + t28 * t20 + pkin(2);
t15 = sin(pkin(6));
t33 = t15 * t17;
t32 = t15 * t20;
t21 = cos(qJ(2));
t31 = t15 * t21;
t30 = cos(pkin(6));
t29 = cos(pkin(10));
t14 = sin(pkin(10));
t27 = t14 * t30;
t26 = t15 * t29;
t25 = t30 * t29;
t23 = t19 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(4) + pkin(8);
t18 = sin(qJ(2));
t9 = t18 * t33 - t30 * t20;
t8 = -t18 * t27 + t29 * t21;
t7 = t29 * t18 + t21 * t27;
t6 = t14 * t21 + t18 * t25;
t5 = t14 * t18 - t21 * t25;
t3 = -t14 * t32 + t8 * t17;
t1 = t6 * t17 + t20 * t26;
t2 = [0, t23 * t8 - t34 * t7, t24 * (t14 * t33 + t8 * t20) - t28 * t3, t3 (-t7 * t16 + t3 * t19) * r_i_i_C(1) + (-t3 * t16 - t7 * t19) * r_i_i_C(2), 0; 0, t23 * t6 - t34 * t5, t24 * (-t17 * t26 + t6 * t20) - t28 * t1, t1 (t1 * t19 - t5 * t16) * r_i_i_C(1) + (-t1 * t16 - t5 * t19) * r_i_i_C(2), 0; 1 (t23 * t18 + t34 * t21) * t15, -t28 * t9 + t24 * (t30 * t17 + t18 * t32) t9 (t16 * t31 + t9 * t19) * r_i_i_C(1) + (-t9 * t16 + t19 * t31) * r_i_i_C(2), 0;];
Ja_transl  = t2;
