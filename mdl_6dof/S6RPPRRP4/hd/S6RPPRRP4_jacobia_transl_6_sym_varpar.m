% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:04
% EndTime: 2019-02-26 20:32:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (129->32), mult. (268->43), div. (0->0), fcn. (351->8), ass. (0->25)
t17 = cos(qJ(4));
t15 = sin(qJ(4));
t30 = pkin(8) + r_i_i_C(2);
t22 = t30 * t15;
t34 = t17 * pkin(4) + pkin(3) + t22;
t33 = pkin(1) + pkin(2);
t14 = sin(qJ(5));
t16 = cos(qJ(5));
t25 = r_i_i_C(3) + qJ(6);
t31 = pkin(5) + r_i_i_C(1);
t32 = t25 * t14 + t31 * t16 + pkin(4);
t29 = cos(qJ(1));
t28 = sin(qJ(1));
t27 = t14 * t17;
t26 = t16 * t17;
t24 = cos(pkin(9));
t23 = sin(pkin(9));
t10 = t29 * t23 - t28 * t24;
t9 = -t28 * t23 - t29 * t24;
t20 = t10 * t26 + t9 * t14;
t19 = t10 * t27 - t9 * t16;
t18 = t32 * t15 - t30 * t17;
t6 = t10 * t14 - t9 * t26;
t5 = -t10 * t16 - t9 * t27;
t1 = [t9 * pkin(7) + t29 * qJ(2) + t34 * t10 + t25 * t19 + t31 * t20 - t33 * t28, t28, 0, t18 * t9, t25 * t6 - t31 * t5, t5; t10 * pkin(7) + t28 * qJ(2) + t25 * t5 + t33 * t29 + t31 * t6 - t34 * t9, -t29, 0, t18 * t10, t19 * t31 - t25 * t20, -t19; 0, 0, -1, -t32 * t17 - t22 (t31 * t14 - t25 * t16) * t15, -t15 * t14;];
Ja_transl  = t1;
