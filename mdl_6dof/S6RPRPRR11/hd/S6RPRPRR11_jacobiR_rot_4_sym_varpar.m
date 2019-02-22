% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:36
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR11_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:36:42
% EndTime: 2019-02-22 10:36:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (69->25), mult. (203->52), div. (0->0), fcn. (284->12), ass. (0->34)
t143 = cos(pkin(6));
t137 = sin(pkin(12));
t147 = cos(qJ(1));
t151 = t147 * t137;
t141 = cos(pkin(12));
t145 = sin(qJ(1));
t152 = t145 * t141;
t131 = t143 * t151 + t152;
t144 = sin(qJ(3));
t146 = cos(qJ(3));
t150 = t147 * t141;
t153 = t145 * t137;
t130 = -t143 * t150 + t153;
t138 = sin(pkin(7));
t142 = cos(pkin(7));
t139 = sin(pkin(6));
t155 = t139 * t147;
t148 = t130 * t142 + t138 * t155;
t123 = -t131 * t146 + t144 * t148;
t157 = t138 * t143;
t156 = t139 * t145;
t154 = t142 * t146;
t149 = t138 * t156;
t122 = -t131 * t144 - t146 * t148;
t140 = cos(pkin(13));
t136 = sin(pkin(13));
t133 = -t143 * t153 + t150;
t132 = -t143 * t152 - t151;
t128 = -t132 * t138 + t142 * t156;
t127 = -t130 * t138 + t142 * t155;
t126 = t146 * t157 + (-t137 * t144 + t141 * t154) * t139;
t125 = t133 * t146 + (t132 * t142 + t149) * t144;
t124 = -t132 * t154 + t133 * t144 - t146 * t149;
t1 = [t123 * t140 + t127 * t136, 0, -t124 * t140, 0, 0, 0; t125 * t140 + t128 * t136, 0, t122 * t140, 0, 0, 0; 0, 0, t126 * t140, 0, 0, 0; -t123 * t136 + t127 * t140, 0, t124 * t136, 0, 0, 0; -t125 * t136 + t128 * t140, 0, -t122 * t136, 0, 0, 0; 0, 0, -t126 * t136, 0, 0, 0; t122, 0, t125, 0, 0, 0; t124, 0, -t123, 0, 0, 0; 0, 0, t144 * t157 + (t141 * t142 * t144 + t137 * t146) * t139, 0, 0, 0;];
JR_rot  = t1;
