% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:16
% EndTime: 2019-02-26 20:56:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (202->34), mult. (167->47), div. (0->0), fcn. (191->10), ass. (0->29)
t18 = cos(qJ(3));
t16 = sin(qJ(3));
t26 = r_i_i_C(2) + qJ(5) + pkin(8);
t20 = t26 * t16;
t17 = cos(qJ(4));
t7 = t17 * pkin(4) + pkin(3);
t32 = t18 * t7 + pkin(2) + t20;
t12 = qJ(4) + pkin(10);
t10 = cos(t12);
t23 = r_i_i_C(3) + qJ(6);
t29 = pkin(5) + r_i_i_C(1);
t8 = sin(t12);
t31 = t29 * t10 + t23 * t8 + t7;
t13 = qJ(1) + pkin(9);
t9 = sin(t13);
t30 = t9 * t8;
t15 = sin(qJ(4));
t28 = pkin(4) * t15;
t27 = t9 * t10;
t11 = cos(t13);
t25 = t11 * t18;
t24 = t15 * t18;
t21 = pkin(7) + t28;
t19 = -t31 * t16 + t26 * t18;
t4 = t10 * t25 + t30;
t3 = t8 * t25 - t27;
t2 = -t11 * t8 + t18 * t27;
t1 = t11 * t10 + t18 * t30;
t5 = [-sin(qJ(1)) * pkin(1) - t29 * t2 + t21 * t11 - t23 * t1 - t32 * t9, 0, t19 * t11, t23 * t4 - t29 * t3 + (-t11 * t24 + t17 * t9) * pkin(4), t11 * t16, t3; cos(qJ(1)) * pkin(1) + t21 * t9 + t29 * t4 + t23 * t3 + t32 * t11, 0, t19 * t9, t23 * t2 - t29 * t1 + (-t11 * t17 - t9 * t24) * pkin(4), t9 * t16, t1; 0, 1, t31 * t18 + t20 (t23 * t10 - t29 * t8 - t28) * t16, -t18, t16 * t8;];
Ja_transl  = t5;
