% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:34:56
% EndTime: 2019-02-22 11:34:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (39->19), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->25)
t139 = sin(qJ(4));
t142 = cos(qJ(2));
t156 = sin(qJ(2));
t157 = cos(qJ(4));
t132 = -t142 * t139 + t156 * t157;
t131 = t156 * t139 + t142 * t157;
t140 = sin(qJ(1));
t127 = t132 * t140;
t138 = sin(qJ(5));
t155 = t127 * t138;
t141 = cos(qJ(5));
t154 = t127 * t141;
t143 = cos(qJ(1));
t130 = t132 * t143;
t153 = t130 * t138;
t152 = t130 * t141;
t151 = t131 * t138;
t150 = t131 * t141;
t128 = t131 * t140;
t146 = -t128 * t138 + t143 * t141;
t144 = t128 * t141 + t143 * t138;
t129 = t131 * t143;
t126 = t129 * t141 - t140 * t138;
t125 = t129 * t138 + t140 * t141;
t1 = [-t144, -t152, 0, t152, -t125, 0; t126, -t154, 0, t154, t146, 0; 0, t150, 0, -t150, -t132 * t138, 0; t127, -t129, 0, t129, 0, 0; -t130, -t128, 0, t128, 0, 0; 0, -t132, 0, t132, 0, 0; t146, -t153, 0, t153, t126, 0; t125, -t155, 0, t155, t144, 0; 0, t151, 0, -t151, t132 * t141, 0;];
JR_rot  = t1;
