% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (192->32), mult. (144->42), div. (0->0), fcn. (157->10), ass. (0->29)
t26 = sin(qJ(6));
t28 = cos(qJ(6));
t52 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
t25 = qJ(2) + pkin(10);
t22 = qJ(4) + t25;
t20 = sin(t22);
t21 = cos(t22);
t38 = pkin(4) + pkin(9) + r_i_i_C(3);
t51 = t20 * qJ(5) + t38 * t21;
t49 = (qJ(5) + t52) * t21;
t40 = pkin(3) * cos(t25) + cos(qJ(2)) * pkin(2);
t48 = pkin(1) + t40 + t51;
t45 = pkin(5) + pkin(8) + qJ(3) + pkin(7);
t27 = sin(qJ(1));
t44 = t27 * t26;
t43 = t27 * t28;
t29 = cos(qJ(1));
t42 = t29 * t26;
t41 = t29 * t28;
t37 = t49 * t27;
t36 = t49 * t29;
t32 = t38 * t20;
t31 = -pkin(3) * sin(t25) - sin(qJ(2)) * pkin(2) - t32;
t30 = t52 * t20 + t51;
t4 = -t20 * t44 + t41;
t3 = t20 * t43 + t42;
t2 = t20 * t42 + t43;
t1 = t20 * t41 - t44;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t48 * t27 + t45 * t29, t31 * t29 + t36, t27, -t29 * t32 + t36, t29 * t20, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t45 * t27 + t48 * t29, t31 * t27 + t37, -t29, -t27 * t32 + t37, t27 * t20, t3 * r_i_i_C(1) + t4 * r_i_i_C(2); 0, t30 + t40, 0, t30, -t21 (-r_i_i_C(1) * t28 + r_i_i_C(2) * t26) * t21;];
Ja_transl  = t5;
