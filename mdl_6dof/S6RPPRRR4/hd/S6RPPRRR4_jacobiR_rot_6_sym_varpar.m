% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (82->17), mult. (118->24), div. (0->0), fcn. (189->8), ass. (0->20)
t121 = cos(qJ(1));
t120 = sin(qJ(1));
t119 = cos(pkin(10));
t118 = sin(pkin(10));
t113 = qJ(5) + qJ(6);
t111 = sin(t113);
t114 = sin(qJ(4));
t106 = t114 * t111;
t112 = cos(t113);
t107 = t114 * t112;
t115 = cos(qJ(4));
t117 = t115 * t111;
t116 = t115 * t112;
t104 = -t120 * t118 - t121 * t119;
t105 = t121 * t118 - t120 * t119;
t101 = t104 * t111 + t105 * t116;
t100 = -t104 * t112 + t105 * t117;
t103 = -t104 * t116 + t105 * t111;
t102 = t104 * t117 + t105 * t112;
t1 = [t101, 0, 0, t104 * t107, t102, t102; t103, 0, 0, t105 * t107, t100, t100; 0, 0, 0, -t116, t106, t106; -t100, 0, 0, -t104 * t106, -t103, -t103; t102, 0, 0, -t105 * t106, t101, t101; 0, 0, 0, t117, t107, t107; t105 * t114, 0, 0, -t104 * t115, 0, 0; -t104 * t114, 0, 0, -t105 * t115, 0, 0; 0, 0, 0, -t114, 0, 0;];
JR_rot  = t1;
