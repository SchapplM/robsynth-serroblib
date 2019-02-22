% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:28:54
% EndTime: 2019-02-22 09:28:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->31), mult. (319->73), div. (0->0), fcn. (454->12), ass. (0->35)
t154 = sin(pkin(6));
t159 = sin(qJ(5));
t171 = t154 * t159;
t162 = cos(qJ(5));
t170 = t154 * t162;
t157 = cos(pkin(6));
t160 = sin(qJ(2));
t169 = t157 * t160;
t163 = cos(qJ(2));
t168 = t157 * t163;
t158 = sin(qJ(6));
t167 = t158 * t162;
t161 = cos(qJ(6));
t166 = t161 * t162;
t153 = sin(pkin(10));
t156 = cos(pkin(10));
t146 = t153 * t160 - t156 * t168;
t147 = t153 * t163 + t156 * t169;
t152 = sin(pkin(11));
t155 = cos(pkin(11));
t165 = -t146 * t155 + t147 * t152;
t148 = t153 * t168 + t156 * t160;
t149 = -t153 * t169 + t156 * t163;
t164 = -t148 * t155 + t149 * t152;
t135 = t146 * t152 + t147 * t155;
t139 = t148 * t152 + t149 * t155;
t145 = (-t152 * t163 + t155 * t160) * t154;
t144 = (t152 * t160 + t155 * t163) * t154;
t141 = t145 * t162 - t157 * t159;
t140 = -t145 * t159 - t157 * t162;
t131 = t139 * t162 - t153 * t171;
t130 = -t139 * t159 - t153 * t170;
t129 = t135 * t162 + t156 * t171;
t128 = -t135 * t159 + t156 * t170;
t1 = [0, -t139 * t158 + t164 * t166, 0, 0, t130 * t161, -t131 * t158 + t161 * t164; 0, -t135 * t158 + t165 * t166, 0, 0, t128 * t161, -t129 * t158 + t161 * t165; 0, t144 * t166 - t145 * t158, 0, 0, t140 * t161, -t141 * t158 + t144 * t161; 0, -t139 * t161 - t164 * t167, 0, 0, -t130 * t158, -t131 * t161 - t158 * t164; 0, -t135 * t161 - t165 * t167, 0, 0, -t128 * t158, -t129 * t161 - t158 * t165; 0, -t144 * t167 - t145 * t161, 0, 0, -t140 * t158, -t141 * t161 - t144 * t158; 0, t164 * t159, 0, 0, t131, 0; 0, t165 * t159, 0, 0, t129, 0; 0, t144 * t159, 0, 0, t141, 0;];
JR_rot  = t1;
