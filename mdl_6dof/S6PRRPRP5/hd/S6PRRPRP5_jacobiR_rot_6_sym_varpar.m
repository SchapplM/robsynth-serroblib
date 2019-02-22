% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:47
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:47:22
% EndTime: 2019-02-22 09:47:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->30), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t138 = sin(pkin(6));
t142 = sin(qJ(3));
t154 = t138 * t142;
t145 = cos(qJ(3));
t153 = t138 * t145;
t140 = cos(pkin(6));
t143 = sin(qJ(2));
t152 = t140 * t143;
t146 = cos(qJ(2));
t151 = t140 * t146;
t141 = sin(qJ(5));
t150 = t141 * t142;
t149 = t141 * t146;
t144 = cos(qJ(5));
t148 = t142 * t144;
t147 = t144 * t146;
t139 = cos(pkin(10));
t137 = sin(pkin(10));
t135 = t140 * t142 + t143 * t153;
t134 = -t140 * t145 + t143 * t154;
t133 = -t137 * t152 + t139 * t146;
t132 = t137 * t151 + t139 * t143;
t131 = t137 * t146 + t139 * t152;
t130 = t137 * t143 - t139 * t151;
t129 = t133 * t145 + t137 * t154;
t128 = t133 * t142 - t137 * t153;
t127 = t131 * t145 - t139 * t154;
t126 = t131 * t142 + t139 * t153;
t1 = [0, -t132 * t150 + t133 * t144, t129 * t141, 0, t128 * t144 - t132 * t141, 0; 0, -t130 * t150 + t131 * t144, t127 * t141, 0, t126 * t144 - t130 * t141, 0; 0 (t142 * t149 + t143 * t144) * t138, t135 * t141, 0, t134 * t144 + t138 * t149, 0; 0, -t132 * t145, -t128, 0, 0, 0; 0, -t130 * t145, -t126, 0, 0, 0; 0, t146 * t153, -t134, 0, 0, 0; 0, t132 * t148 + t133 * t141, -t129 * t144, 0, t128 * t141 + t132 * t144, 0; 0, t130 * t148 + t131 * t141, -t127 * t144, 0, t126 * t141 + t130 * t144, 0; 0 (t141 * t143 - t142 * t147) * t138, -t135 * t144, 0, t134 * t141 - t138 * t147, 0;];
JR_rot  = t1;
