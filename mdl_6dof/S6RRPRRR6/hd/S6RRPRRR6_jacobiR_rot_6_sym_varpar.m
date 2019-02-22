% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:42:56
% EndTime: 2019-02-22 11:42:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (136->24), mult. (170->24), div. (0->0), fcn. (256->8), ass. (0->26)
t150 = qJ(4) + qJ(5);
t139 = sin(t150);
t143 = cos(qJ(2));
t148 = cos(t150);
t155 = sin(qJ(2));
t134 = -t143 * t139 + t155 * t148;
t133 = t155 * t139 + t143 * t148;
t156 = cos(qJ(1));
t141 = sin(qJ(1));
t128 = t134 * t141;
t140 = sin(qJ(6));
t154 = t128 * t140;
t142 = cos(qJ(6));
t126 = t128 * t142;
t131 = t134 * t156;
t153 = t131 * t140;
t127 = t131 * t142;
t152 = t133 * t140;
t132 = t133 * t142;
t129 = t133 * t141;
t145 = t129 * t140 - t156 * t142;
t144 = -t129 * t142 - t156 * t140;
t130 = t133 * t156;
t125 = t130 * t142 - t141 * t140;
t124 = -t130 * t140 - t141 * t142;
t1 = [t144, -t127, 0, t127, t127, t124; t125, -t126, 0, t126, t126, -t145; 0, t132, 0, -t132, -t132, -t134 * t140; t145, t153, 0, -t153, -t153, -t125; t124, t154, 0, -t154, -t154, t144; 0, -t152, 0, t152, t152, -t134 * t142; t128, -t130, 0, t130, t130, 0; -t131, -t129, 0, t129, t129, 0; 0, -t134, 0, t134, t134, 0;];
JR_rot  = t1;
