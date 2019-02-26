% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP7_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:55
% EndTime: 2019-02-26 22:28:55
% DurationCPUTime: 0.09s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t104 = sin(pkin(6));
t108 = sin(qJ(1));
t117 = t108 * t104;
t107 = sin(qJ(2));
t116 = t108 * t107;
t110 = cos(qJ(2));
t115 = t108 * t110;
t111 = cos(qJ(1));
t114 = t111 * t104;
t113 = t111 * t107;
t112 = t111 * t110;
t109 = cos(qJ(3));
t106 = sin(qJ(3));
t105 = cos(pkin(6));
t1 = [0, t117, t105 * t115 + t113 (-t105 * t116 + t112) * t106 - t109 * t117, 0, 0; 0, -t114, -t105 * t112 + t116 (t105 * t113 + t115) * t106 + t109 * t114, 0, 0; 1, t105, -t104 * t110, t104 * t107 * t106 - t105 * t109, 0, 0;];
Jg_rot  = t1;
