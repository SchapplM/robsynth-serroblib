% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR8_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:02
% EndTime: 2019-02-22 11:17:02
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t62 = sin(qJ(2));
t63 = sin(qJ(1));
t70 = t63 * t62;
t60 = sin(pkin(10));
t64 = cos(qJ(2));
t69 = t64 * t60;
t61 = cos(pkin(10));
t68 = t64 * t61;
t65 = cos(qJ(1));
t67 = t65 * t62;
t66 = t65 * t64;
t1 = [t65 * t60 - t63 * t68, -t61 * t67, 0, 0, 0, 0; t63 * t60 + t61 * t66, -t61 * t70, 0, 0, 0, 0; 0, t68, 0, 0, 0, 0; -t70, t66, 0, 0, 0, 0; t67, t63 * t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; -t65 * t61 - t63 * t69, -t60 * t67, 0, 0, 0, 0; t60 * t66 - t63 * t61, -t60 * t70, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0;];
JR_rot  = t1;
