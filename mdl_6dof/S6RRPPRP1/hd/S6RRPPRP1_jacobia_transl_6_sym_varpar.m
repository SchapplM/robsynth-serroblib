% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:50
% EndTime: 2019-02-26 21:24:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (188->31), mult. (160->40), div. (0->0), fcn. (186->10), ass. (0->25)
t15 = qJ(2) + pkin(9);
t10 = sin(t15);
t29 = r_i_i_C(2) + pkin(8) + qJ(4);
t34 = cos(qJ(2)) * pkin(2) + t29 * t10;
t12 = cos(t15);
t7 = cos(pkin(10)) * pkin(4) + pkin(3);
t33 = t12 * t7 + pkin(1) + t34;
t14 = pkin(10) + qJ(5);
t11 = cos(t14);
t26 = r_i_i_C(3) + qJ(6);
t31 = pkin(5) + r_i_i_C(1);
t9 = sin(t14);
t32 = t31 * t11 + t26 * t9 + t7;
t20 = sin(qJ(1));
t30 = t20 * t9;
t21 = cos(qJ(1));
t28 = t12 * t21;
t27 = t20 * t11;
t23 = pkin(4) * sin(pkin(10)) + qJ(3) + pkin(7);
t22 = -sin(qJ(2)) * pkin(2) + t29 * t12 - t32 * t10;
t4 = t11 * t28 + t30;
t3 = t9 * t28 - t27;
t2 = t12 * t27 - t21 * t9;
t1 = t11 * t21 + t12 * t30;
t5 = [-t26 * t1 - t31 * t2 - t33 * t20 + t23 * t21, t22 * t21, t20, t21 * t10, t26 * t4 - t31 * t3, t3; t23 * t20 + t33 * t21 + t26 * t3 + t31 * t4, t22 * t20, -t21, t20 * t10, -t31 * t1 + t26 * t2, t1; 0, t32 * t12 + t34, 0, -t12 (t26 * t11 - t31 * t9) * t10, t10 * t9;];
Ja_transl  = t5;
