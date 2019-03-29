% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobia_transl_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:50
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->25), mult. (117->40), div. (0->0), fcn. (129->10), ass. (0->28)
t21 = qJ(3) + qJ(4);
t17 = sin(t21);
t16 = t17 * r_i_i_C(3);
t40 = cos(qJ(3)) * pkin(2);
t42 = t16 + t40;
t19 = cos(t21);
t23 = sin(qJ(5));
t41 = r_i_i_C(2) * t17 * t23 + r_i_i_C(3) * t19;
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t38 = t41 * t18;
t37 = t19 * t23;
t25 = cos(qJ(5));
t36 = t19 * t25;
t20 = cos(t22);
t35 = t20 * t23;
t34 = t20 * t25;
t33 = t41 * t20;
t32 = r_i_i_C(1) * t17 * t25;
t7 = t18 * t25 - t19 * t35;
t8 = t18 * t23 + t19 * t34;
t30 = t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t42 * t20;
t29 = r_i_i_C(1) * t36 - r_i_i_C(2) * t37 + t16;
t28 = -sin(qJ(3)) * pkin(2) - t32;
t5 = t18 * t37 + t34;
t6 = -t18 * t36 + t35;
t27 = t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t18;
t1 = [-sin(qJ(1)) * pkin(1) + t27, t27, t28 * t20 + t33, -t20 * t32 + t33, r_i_i_C(1) * t7 - r_i_i_C(2) * t8; cos(qJ(1)) * pkin(1) + t30, t30, t28 * t18 + t38, -t18 * t32 + t38, -r_i_i_C(1) * t5 + r_i_i_C(2) * t6; 0, 0, t29 + t40, t29 (-r_i_i_C(1) * t23 - r_i_i_C(2) * t25) * t17;];
Ja_transl  = t1;
