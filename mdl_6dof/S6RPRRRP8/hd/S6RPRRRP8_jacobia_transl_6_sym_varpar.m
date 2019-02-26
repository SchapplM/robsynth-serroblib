% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:04
% EndTime: 2019-02-26 21:12:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (160->33), mult. (200->41), div. (0->0), fcn. (223->8), ass. (0->29)
t17 = qJ(3) + qJ(4);
t15 = sin(t17);
t43 = pkin(9) + r_i_i_C(2);
t45 = t15 * t43;
t16 = cos(t17);
t44 = t43 * t16;
t42 = pkin(5) + r_i_i_C(1);
t32 = r_i_i_C(3) + qJ(6);
t40 = cos(qJ(3)) * pkin(3);
t39 = sin(qJ(3)) * pkin(3);
t38 = pkin(1) + pkin(8) + pkin(7);
t18 = sin(qJ(5));
t20 = sin(qJ(1));
t36 = t20 * t18;
t21 = cos(qJ(5));
t35 = t20 * t21;
t23 = cos(qJ(1));
t34 = t23 * t18;
t33 = t23 * t21;
t29 = t20 * t45 + (pkin(4) * t20 + t32 * t36 + t35 * t42) * t16;
t28 = -t32 * t18 - t42 * t21 - pkin(4);
t27 = t15 * pkin(4) + qJ(2) + t39 - t44;
t26 = t28 * t15 + t44;
t25 = t28 * t16 - t45;
t4 = t15 * t33 - t36;
t3 = t15 * t34 + t35;
t2 = t15 * t35 + t34;
t1 = t15 * t36 - t33;
t5 = [-t38 * t20 + t27 * t23 + t32 * t3 + t42 * t4, t20, t20 * t40 + t29, t29, -t42 * t1 + t32 * t2, t1; t32 * t1 + t42 * t2 + t27 * t20 + t38 * t23, -t23 (t25 - t40) * t23, t25 * t23, t42 * t3 - t32 * t4, -t3; 0, 0, t26 - t39, t26 (-t42 * t18 + t32 * t21) * t16, t16 * t18;];
Ja_transl  = t5;
