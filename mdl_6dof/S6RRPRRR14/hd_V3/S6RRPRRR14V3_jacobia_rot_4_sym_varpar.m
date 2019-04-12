% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14V3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
t39 = cos(qJ(2));
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t45 = t37 * t36;
t30 = atan2(-t45, -t39);
t28 = sin(t30);
t29 = cos(t30);
t21 = -t28 * t45 - t29 * t39;
t20 = 0.1e1 / t21 ^ 2;
t40 = cos(qJ(1));
t50 = t20 * t40 ^ 2;
t35 = sin(qJ(4));
t38 = cos(qJ(4));
t42 = t40 * t38;
t27 = t37 * t35 + t39 * t42;
t25 = 0.1e1 / t27 ^ 2;
t43 = t40 * t35;
t26 = -t37 * t38 + t39 * t43;
t49 = t25 * t26;
t32 = t36 ^ 2;
t48 = t32 / t39 ^ 2;
t47 = t36 * t40;
t31 = 0.1e1 / (t37 ^ 2 * t48 + 0.1e1);
t46 = t37 * t31;
t44 = t37 * t39;
t41 = t26 ^ 2 * t25 + 0.1e1;
t33 = 0.1e1 / t39;
t24 = 0.1e1 / t27;
t23 = (0.1e1 + t48) * t46;
t22 = 0.1e1 / t41;
t19 = 0.1e1 / t21;
t18 = 0.1e1 / (t32 * t50 + 0.1e1);
t1 = [t33 * t31 * t47, t23, 0, 0, 0, 0; (-t19 * t45 - (-t29 * t32 * t33 * t46 + (t31 - 0.1e1) * t36 * t28) * t36 * t50) * t18 (t39 * t19 - (-t28 * t44 + t29 * t36 + (t28 * t39 - t29 * t45) * t23) * t36 * t20) * t40 * t18, 0, 0, 0, 0; ((-t35 * t44 - t42) * t24 - (-t38 * t44 + t43) * t49) * t22 (-t24 * t35 + t38 * t49) * t22 * t47, 0, t41 * t22, 0, 0;];
Ja_rot  = t1;
