% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:05
% EndTime: 2019-02-26 20:32:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->15), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
t115 = cos(qJ(1));
t114 = sin(qJ(1));
t113 = cos(pkin(9));
t112 = sin(pkin(9));
t102 = sin(qJ(5));
t103 = sin(qJ(4));
t111 = t103 * t102;
t104 = cos(qJ(5));
t110 = t103 * t104;
t105 = cos(qJ(4));
t109 = t105 * t102;
t108 = t105 * t104;
t97 = -t114 * t112 - t115 * t113;
t98 = t115 * t112 - t114 * t113;
t107 = t97 * t102 + t98 * t108;
t106 = -t97 * t104 + t98 * t109;
t96 = t98 * t102 - t97 * t108;
t95 = -t98 * t104 - t97 * t109;
t1 = [t107, 0, 0, t97 * t110, -t95, 0; t96, 0, 0, t98 * t110, t106, 0; 0, 0, 0, -t108, t111, 0; t98 * t103, 0, 0, -t97 * t105, 0, 0; -t97 * t103, 0, 0, -t98 * t105, 0, 0; 0, 0, 0, -t103, 0, 0; t106, 0, 0, t97 * t111, t96, 0; t95, 0, 0, t98 * t111, -t107, 0; 0, 0, 0, -t109, -t110, 0;];
JR_rot  = t1;
