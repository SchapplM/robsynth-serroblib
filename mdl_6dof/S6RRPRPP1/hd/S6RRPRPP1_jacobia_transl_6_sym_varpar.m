% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:51
% EndTime: 2019-02-26 21:34:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (196->35), mult. (172->47), div. (0->0), fcn. (198->10), ass. (0->29)
t15 = qJ(2) + pkin(9);
t10 = sin(t15);
t32 = r_i_i_C(2) + qJ(5) + pkin(8);
t37 = cos(qJ(2)) * pkin(2) + t32 * t10;
t12 = cos(t15);
t21 = cos(qJ(4));
t7 = t21 * pkin(4) + pkin(3);
t36 = t12 * t7 + pkin(1) + t37;
t14 = qJ(4) + pkin(10);
t11 = cos(t14);
t27 = r_i_i_C(3) + qJ(6);
t34 = pkin(5) + r_i_i_C(1);
t9 = sin(t14);
t35 = t34 * t11 + t27 * t9 + t7;
t18 = sin(qJ(4));
t33 = pkin(4) * t18;
t20 = sin(qJ(1));
t31 = t12 * t20;
t22 = cos(qJ(1));
t30 = t12 * t22;
t29 = t20 * t11;
t28 = t22 * t11;
t24 = qJ(3) + pkin(7) + t33;
t23 = -sin(qJ(2)) * pkin(2) + t32 * t12 - t35 * t10;
t4 = t12 * t28 + t20 * t9;
t3 = t9 * t30 - t29;
t2 = t12 * t29 - t22 * t9;
t1 = t9 * t31 + t28;
t5 = [-t27 * t1 - t34 * t2 - t36 * t20 + t24 * t22, t23 * t22, t20, t27 * t4 - t34 * t3 + (-t18 * t30 + t20 * t21) * pkin(4), t22 * t10, t3; t24 * t20 + t36 * t22 + t27 * t3 + t34 * t4, t23 * t20, -t22, t27 * t2 - t34 * t1 + (-t18 * t31 - t21 * t22) * pkin(4), t20 * t10, t1; 0, t35 * t12 + t37, 0 (t27 * t11 - t34 * t9 - t33) * t10, -t12, t10 * t9;];
Ja_transl  = t5;
