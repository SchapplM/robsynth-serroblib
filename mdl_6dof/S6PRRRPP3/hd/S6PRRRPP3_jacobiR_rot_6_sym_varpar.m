% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:53
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:53:31
% EndTime: 2019-02-22 09:53:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t133 = sin(pkin(6));
t137 = sin(qJ(3));
t149 = t133 * t137;
t140 = cos(qJ(3));
t148 = t133 * t140;
t141 = cos(qJ(2));
t147 = t133 * t141;
t135 = cos(pkin(6));
t138 = sin(qJ(2));
t146 = t135 * t138;
t145 = t135 * t141;
t136 = sin(qJ(4));
t144 = t136 * t140;
t139 = cos(qJ(4));
t143 = t139 * t140;
t142 = t140 * t141;
t134 = cos(pkin(10));
t132 = sin(pkin(10));
t130 = t135 * t137 + t138 * t148;
t129 = t135 * t140 - t138 * t149;
t128 = -t132 * t146 + t134 * t141;
t127 = t132 * t145 + t134 * t138;
t126 = t132 * t141 + t134 * t146;
t125 = t132 * t138 - t134 * t145;
t124 = t128 * t140 + t132 * t149;
t123 = -t128 * t137 + t132 * t148;
t122 = t126 * t140 - t134 * t149;
t121 = -t126 * t137 - t134 * t148;
t1 = [0, -t127 * t137, t124, 0, 0, 0; 0, -t125 * t137, t122, 0, 0, 0; 0, t137 * t147, t130, 0, 0, 0; 0, -t127 * t144 - t128 * t139, t123 * t136, t124 * t139 + t127 * t136, 0, 0; 0, -t125 * t144 - t126 * t139, t121 * t136, t122 * t139 + t125 * t136, 0, 0; 0 (t136 * t142 - t138 * t139) * t133, t129 * t136, t130 * t139 - t136 * t147, 0, 0; 0, -t127 * t143 + t128 * t136, t123 * t139, -t124 * t136 + t127 * t139, 0, 0; 0, -t125 * t143 + t126 * t136, t121 * t139, -t122 * t136 + t125 * t139, 0, 0; 0 (t136 * t138 + t139 * t142) * t133, t129 * t139, -t130 * t136 - t139 * t147, 0, 0;];
JR_rot  = t1;
