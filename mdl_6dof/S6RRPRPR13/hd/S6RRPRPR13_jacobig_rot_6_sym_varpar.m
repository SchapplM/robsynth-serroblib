% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (8->8), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t103 = sin(pkin(6));
t107 = sin(qJ(1));
t116 = t107 * t103;
t106 = sin(qJ(2));
t115 = t107 * t106;
t109 = cos(qJ(2));
t114 = t107 * t109;
t110 = cos(qJ(1));
t113 = t110 * t103;
t112 = t110 * t106;
t111 = t110 * t109;
t108 = cos(qJ(4));
t105 = sin(qJ(4));
t104 = cos(pkin(6));
t1 = [0, t116, 0, -t104 * t115 + t111, 0, t105 * t116 - (t104 * t114 + t112) * t108; 0, -t113, 0, t104 * t112 + t114, 0, -t105 * t113 - (-t104 * t111 + t115) * t108; 1, t104, 0, t103 * t106, 0, t103 * t109 * t108 + t104 * t105;];
Jg_rot  = t1;
