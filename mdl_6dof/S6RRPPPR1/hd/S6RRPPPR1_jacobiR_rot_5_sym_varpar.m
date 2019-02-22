% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:06
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:06:18
% EndTime: 2019-02-22 11:06:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t73 = sin(pkin(10));
t75 = sin(qJ(1));
t80 = t75 * t73;
t74 = cos(pkin(10));
t79 = t75 * t74;
t76 = cos(qJ(1));
t78 = t76 * t73;
t77 = t76 * t74;
t72 = qJ(2) + pkin(9);
t71 = cos(t72);
t70 = sin(t72);
t1 = [-t71 * t79 + t78, -t70 * t77, 0, 0, 0, 0; t71 * t77 + t80, -t70 * t79, 0, 0, 0, 0; 0, t71 * t74, 0, 0, 0, 0; -t75 * t70, t76 * t71, 0, 0, 0, 0; t76 * t70, t75 * t71, 0, 0, 0, 0; 0, t70, 0, 0, 0, 0; -t71 * t80 - t77, -t70 * t78, 0, 0, 0, 0; t71 * t78 - t79, -t70 * t80, 0, 0, 0, 0; 0, t71 * t73, 0, 0, 0, 0;];
JR_rot  = t1;
