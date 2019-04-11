% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10V2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:36
% EndTime: 2019-04-11 14:56:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->25), mult. (119->37), div. (0->0), fcn. (127->8), ass. (0->28)
t19 = qJ(2) + qJ(3);
t16 = sin(t19);
t17 = cos(t19);
t20 = sin(qJ(4));
t38 = r_i_i_C(2) * t20;
t44 = pkin(5) + r_i_i_C(3);
t45 = t16 * t38 + t17 * t44;
t42 = t17 * pkin(3) + t44 * t16;
t18 = cos(qJ(2)) * pkin(2);
t41 = pkin(1) + t18 + t42;
t23 = cos(qJ(4));
t39 = r_i_i_C(1) * t23;
t24 = cos(qJ(1));
t35 = t20 * t24;
t22 = sin(qJ(1));
t34 = t22 * t20;
t33 = t22 * t23;
t32 = t23 * t24;
t31 = t45 * t22;
t29 = t45 * t24;
t27 = (-pkin(3) - t39) * t16;
t26 = (-t38 + t39) * t17 + t42;
t25 = -sin(qJ(2)) * pkin(2) + t27;
t4 = t17 * t32 + t34;
t3 = -t17 * t35 + t33;
t2 = -t17 * t33 + t35;
t1 = t17 * t34 + t32;
t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t41 * t22, t25 * t24 + t29, t24 * t27 + t29, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t41 * t24, t25 * t22 + t31, t22 * t27 + t31, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t26, t26 (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
Ja_transl  = t5;
