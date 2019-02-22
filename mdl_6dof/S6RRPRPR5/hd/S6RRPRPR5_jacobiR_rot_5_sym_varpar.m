% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:25:09
% EndTime: 2019-02-22 11:25:09
% DurationCPUTime: 0.14s
% Computational Cost: add. (114->27), mult. (317->61), div. (0->0), fcn. (452->12), ass. (0->32)
t146 = sin(pkin(12));
t155 = cos(qJ(4));
t165 = t146 * t155;
t148 = sin(pkin(6));
t154 = sin(qJ(1));
t164 = t148 * t154;
t157 = cos(qJ(1));
t163 = t148 * t157;
t149 = cos(pkin(12));
t162 = t149 * t155;
t151 = cos(pkin(6));
t147 = sin(pkin(11));
t150 = cos(pkin(11));
t153 = sin(qJ(2));
t156 = cos(qJ(2));
t159 = t156 * t147 + t153 * t150;
t141 = t159 * t151;
t142 = t153 * t147 - t156 * t150;
t161 = t157 * t141 - t154 * t142;
t160 = -t154 * t141 - t157 * t142;
t152 = sin(qJ(4));
t127 = -t152 * t161 - t155 * t163;
t128 = t152 * t163 - t155 * t161;
t158 = t142 * t151;
t140 = t159 * t148;
t139 = t142 * t148;
t137 = -t140 * t152 + t151 * t155;
t135 = t154 * t158 - t157 * t159;
t132 = -t154 * t159 - t157 * t158;
t130 = t152 * t164 + t155 * t160;
t129 = t152 * t160 - t155 * t164;
t1 = [t128 * t149 + t132 * t146, t135 * t162 + t146 * t160, 0, -t129 * t149, 0, 0; t130 * t149 - t135 * t146, t132 * t162 + t146 * t161, 0, t127 * t149, 0, 0; 0, -t139 * t162 + t140 * t146, 0, t137 * t149, 0, 0; -t128 * t146 + t132 * t149, -t135 * t165 + t149 * t160, 0, t129 * t146, 0, 0; -t130 * t146 - t135 * t149, -t132 * t165 + t149 * t161, 0, -t127 * t146, 0, 0; 0, t139 * t165 + t140 * t149, 0, -t137 * t146, 0, 0; t127, t135 * t152, 0, t130, 0, 0; t129, t132 * t152, 0, -t128, 0, 0; 0, -t139 * t152, 0, t140 * t155 + t151 * t152, 0, 0;];
JR_rot  = t1;
